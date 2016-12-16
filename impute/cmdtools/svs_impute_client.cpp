/* Copyright 2016 Golden Helix, Inc. */
#include <QBuffer>
#include <QDataStream>
#include <QCoreApplication>
#include <QDate>
#include <QDir>
#include <QFile>
#include <QStringList>
#include <QTextStream>
#include <QTimer>
#include <QLocalSocket>
#include <QThread>

#include "impute/cmdtools/imputeopts.h"

// ---------------
// Control Channel
// ---------------
class ControlStream
{
public:
  ControlStream(QFileDevice* control) : _control(control) {}
  void sendMessage(QString msg, QStringList args)
  {
    _control->write(toC(msg));
    foreach(QString arg, args){
      _control->write("\t");
      _control->write(toC(arg));
    }
    _control->write("\n");
    // Important, all IO is buffered by default. Flush actually pushes
    // buffered data down through OS channels and across process
    // boundaries!
    _control->flush();
  }

  // Convenience methods
  void sendMessage(QString msg) { sendMessage(msg, QStringList()); }
  void sendMessage(QString msg, QString arg1) { sendMessage(msg, QStringList() << arg1); }
  void sendMessage(QString msg, QString arg1, QString arg2) { sendMessage(msg, QStringList() << arg1 << arg2); }
  void sendMessage(QString msg, QString arg1, QString arg2, QString arg3) { sendMessage(msg, QStringList() << arg1 << arg2 << arg3); }

  void sendError(QString errMsg){
    sendMessage("ERROR", errMsg);
  }
  void sendDebug(QString errMsg){
    sendMessage("DEBUG_OUTPUT", errMsg);
  }

private:
  QFileDevice* _control;
};

// Wrap a QLocalServer, but block on readData using waitForReadyRead
// until buffer has size or connection is lost.
class SyncLockSocketWrapper : public QIODevice
{
  Q_OBJECT;

public:
  SyncLockSocketWrapper(QObject* parent = 0) : QIODevice(parent), _dev(0), _bytesRead(0) {}

  // Don't use this class until you set the backend device
  void setDevice(QLocalSocket* dev){
    _dev = dev;
    // Pass through signals
    connect(_dev, SIGNAL(aboutToClose()), this, SIGNAL(aboutToClose()));
    connect(_dev, SIGNAL(bytesWritten(qint64)), this, SIGNAL(bytesWritten(qint64)));
    connect(_dev, SIGNAL(readChannelFinished()), this, SIGNAL(readChannelFinished()));
    connect(_dev, SIGNAL(readyRead()), this, SIGNAL(readyRead()));
  }

  void clearBuffers(){
    QIODevice::close();
    QIODevice::open(QIODevice::ReadWrite);
  }

  // QIODevice does buffering, so we need to account for our own
  // instance buffer here
  qint64 bytesAvailable() const {
    qint64 available = QIODevice::bytesAvailable();
    available += _dev->bytesAvailable();
    return available;
  }

  // Pass through
  bool isSequential() const { return _dev->isSequential(); }
  qint64 bytesToWrite() const { return _dev->bytesToWrite(); }
  bool canReadLine() const { return _dev->canReadLine(); }
  void close() const { return _dev->close(); }
  bool waitForBytesWritten(int msecs = 30000) { return _dev->waitForBytesWritten(msecs); }
  bool waitForReadyRead(int msecs = 30000) { return _dev->waitForReadyRead(msecs); }
protected:
  qint64 readData(char* data, qint64 size)
  {
    if(!size)
      _dev->read(data, size);
    qint64 sizeLeft = size;
    char* dataPtr = data;
    while (sizeLeft > 0){
      qint64 avail = _dev->bytesAvailable();
      if(avail == 0){
        if(!_dev->waitForReadyRead(10)) {
          if (_dev->state() != QLocalSocket::ConnectedState)
            return (dataPtr-data);

          // Socket closing is a queued connection, run an event loop
          // to see if there is more data or the socket has reached
          // its end and will close itself.
          qDebug("Waiting for end of stream...");
          QEventLoop wait;
          connect(_dev, SIGNAL(stateChanged(QLocalSocket::LocalSocketState)), &wait, SLOT(quit()));
          connect(_dev, SIGNAL(readyRead()), &wait, SLOT(quit()));
          wait.exec();
        }
        continue;
      }

      qint64 ret = _dev->read(dataPtr, sizeLeft);
      Q_ASSERT(ret > 0 && ret <= size && ret <= avail);
      sizeLeft -= ret;
      dataPtr = dataPtr + (int)ret;
      _bytesRead += ret;
      //qDebug("read %lld total %lld of %lld avail (%lld left)", _bytesRead, size, avail, sizeLeft);
    }
    return size;
  }
  // Pass through
  qint64 writeData(const char* data, qint64 size) { return _dev->write(data, size); }
private:
  QLocalSocket* _dev;
  qint64 _bytesRead;
};

// ---------------
// StreamDataParser
// ---------------
class StreamDataParser
{
public:
  StreamDataParser(QString pipeName, ControlStream& control, QString inputToken);
  bool readSampleNames();
  QList<QByteArray> sampleNames() const { return _sampleNames; }

  // Contains a full variant's data for all sample (QVector instances
  // will match _nSamples)
  struct Record{
    QByteArray chrom;
    int pos;
    QByteArray id;
    QList<QByteArray> alleles;

    QVector<int> var1;
    QVector<int> var2;
    QVector<bool> phased;
  };

  bool hasNextRec();
  Record takeNextRec() { return _recordBuffer.takeFirst(); }

private:
  int _nSamples;
  QString _inputToken;
  QList<QByteArray> _sampleNames;
  QList<Record> _recordBuffer;

  QString _pipeName;
  ControlStream& _control;
  QLocalSocket _socket;
  SyncLockSocketWrapper _syncSocket;
  QDataStream _input;

};

StreamDataParser::StreamDataParser(QString pipeName, ControlStream& control, QString inputToken)
  : _pipeName(pipeName), _control(control), _inputToken(inputToken)
{
  _input.setDevice(&_syncSocket);
  _input.setByteOrder(QDataStream::LittleEndian);
  _syncSocket.setDevice(&_socket);
}

bool StreamDataParser::readSampleNames()
{
  _control.sendMessage("READ_SAMPLES", _inputToken);
  _socket.connectToServer(_pipeName, QIODevice::ReadOnly);
  if(!_socket.waitForConnected(5000)){
    _control.sendError("Unable to connect to pip named: " + _pipeName);
    exit(1);
  }
  _syncSocket.clearBuffers();

  _input >> _nSamples;
  if(_nSamples == 0)
    return false;

  for(int i=0; i<_nSamples; i++){
    QByteArray sample;
    _input >> sample;
    _sampleNames << sample;
  }
  _control.sendDebug(QString("Received %1 samples").arg(_nSamples));
  return true;
}

bool StreamDataParser::hasNextRec()
{
  if(_recordBuffer.size() > 0)
    return true;

  // Fill up record buffer
  _control.sendMessage("READ_RECORDS", _inputToken);

  _socket.connectToServer(_pipeName, QIODevice::ReadOnly);
  if(!_socket.waitForConnected(5000)){
    _control.sendError("Unable to connect to pip named: " + _pipeName);
    exit(1);
  }
  _syncSocket.clearBuffers();

  int numRecords;
  _input >> numRecords;

  // Don't bother asking for more data
  if(numRecords == 0) {
    return false;
  }

  // Should be ready to read numRecords into our buffer
  for(int r=0; r<numRecords; r++){
    qDebug("reading %d of %d", r, numRecords);
    Record record;
    _input >> record.chrom;
    _input >> record.pos;
    _input >> record.id;
    _input >> record.alleles;
    _input >> record.var1;
    _input >> record.var2;
    _input >> record.phased;

    _recordBuffer << record;
  }
  return true;
}

// -----------------
// StreamRefDataReader
// -----------------
class StreamRefDataReader : public RefDataReader
{
public:
  StreamRefDataReader(QString pipeName, ControlStream& control, QString inputToken);

  bool canAdvanceWindow() const { return _nextRecordExists; }
  bool hasNextRec() const { return _nextRecordExists; }
  BitSetRefGT nextRec() const;

  void advanceRec();

  // void addNewDataToNewWindow(int windowSize);   // Use default implementation here.
  // bool lastWindowOnChrom() const;               // Use default implementation here.

private:
  void assembleNextRefRec();

  int _nSamples;
  BitSetRefGT _nextRecord;
  bool _nextRecordExists;
  StreamDataParser _parser;
};

StreamRefDataReader::StreamRefDataReader(QString pipeName, ControlStream& control, QString inputToken)
  : _parser(pipeName, control, inputToken)
{
  if (!_parser.readSampleNames()) {
    control.sendError("Expected sample names");
    exit(1);
  }

  foreach (QByteArray sampName, _parser.sampleNames())
    _samples.setSamp(SampleNames::getIndex(sampName));

  _nSamples = _samples.nSamples();

  // Cache the first record, assuming it exists. Else mark that we
  // don't have a "next" record. Also set up the "samples" object.
  _nextRecordExists = _parser.hasNextRec();
  if(_nextRecordExists)
    assembleNextRefRec();
}

BitSetRefGT StreamRefDataReader::nextRec() const
{
  return _nextRecord;
}

void StreamRefDataReader::advanceRec()
{
  // Read the record AFTER the "next" record. (The "next" record is
  // the one already in memory.) Mark if it doesn't read, which should
  // usually happen only if we are at the EOF.

  _nextRecordExists = _parser.hasNextRec();

  if (_nextRecordExists)
    assembleNextRefRec();
}

void StreamRefDataReader::assembleNextRefRec()
{
  StreamDataParser::Record record = _parser.takeNextRec();
  BitSetRefGT newrr(_samples);

  newrr.setIdInfo(ChromeIds::getIndex(record.chrom), record.pos, record.id);

  foreach (QByteArray al, record.alleles)
    newrr.addAllele(al);

  newrr.storePhasedAlleles(record.var1, record.var2);

  _nextRecord = newrr;
}

// -----------------
// StreamTargetDataReader
// -----------------
class StreamTargetDataReader : public TargDataReader
{
public:
  StreamTargetDataReader(QString pipeName, ControlStream& control, QString inputToken);

  bool canAdvanceWindow() const { return _nextRecordExists; }
  bool hasNextRec() const { return _nextRecordExists; }
  BitSetGT nextRec() const;

  void advanceRec();

  // void addNewDataToNewWindow(int windowSize);   // Use default implementation here.
  // bool lastWindowOnChrom() const;               // Use default implementation here.

private:
  void assembleNextTargetRec();

  BitSetGT _nextRecord;
  bool _nextRecordExists;
  StreamDataParser _parser;
};

StreamTargetDataReader::StreamTargetDataReader(QString pipeName, ControlStream& control,
                                               QString inputToken)
  : _parser(pipeName, control, inputToken)
{
  if (!_parser.readSampleNames()) {
    control.sendError("Expected sample names");
    exit(1);
  }
  foreach (QByteArray sampName, _parser.sampleNames())
    _samples.setSamp(SampleNames::getIndex(sampName));

  // Cache the first record, assuming it exists. Else mark that we
  // don't have a "next" record. Also set up the "samples" object.
  _nextRecordExists = _parser.hasNextRec();
  if(_nextRecordExists)
    assembleNextTargetRec();
}

BitSetGT StreamTargetDataReader::nextRec() const
{
  return _nextRecord;
}

void StreamTargetDataReader::advanceRec()
{
  // Read the record AFTER the "next" record. (The "next" record is
  // the one already in memory.) Mark if it doesn't read, which should
  // usually happen only if we are at the EOF.
  _nextRecordExists = _parser.hasNextRec();

  if (_nextRecordExists)
    assembleNextTargetRec();
}

void StreamTargetDataReader::assembleNextTargetRec()
{
  StreamDataParser::Record record = _parser.takeNextRec();
  BitSetGT newtr(_samples);

  newtr.setIdInfo(ChromeIds::getIndex(record.chrom), record.pos, record.id);

  foreach (QByteArray al, record.alleles)
    newtr.addAllele(al);

  newtr.storeAlleles(record.var1, record.var2, record.phased);

  _nextRecord = newtr;
}

// --------------
// StreamDataWriter
// --------------
class StreamDataWriter : public ImputeDataWriter
{
public:
  StreamDataWriter(QString pipeName, ControlStream& control, const Samples& samples);
  void writeHeader();
  void writeEOF() {}
protected:
  // These are all overrides of ImputeDataWriter virtual functions
  void initializeWindowBuffering(const int initSize);
  void appendPhasedVariantData();
  void finishAndWriteRec();
  void finalizeForWindow();

private:
  void outputMarker();
  void outputInfo();
  void outputFormat();

  QString _pipeName;
  ControlStream& _control;
  QDataStream _out;
  QLocalSocket _socket;
};

StreamDataWriter::StreamDataWriter(QString pipeName, ControlStream& control, const Samples& samples)
  : ImputeDataWriter(samples), _pipeName(pipeName), _control(control)
{
}

void StreamDataWriter::writeHeader()
{
  // NO-OP
}

void StreamDataWriter::initializeWindowBuffering(const int initSize)
{
  _control.sendDebug(QString("Init window %1").arg(initSize));
  // connect;
  _control.sendMessage("WRITE_RECORDS");

  _socket.connectToServer(_pipeName, QIODevice::WriteOnly);
  if(!_socket.waitForConnected(5000)){
    _control.sendError("Unable to connect to pip named: " + _pipeName);
    exit(1);
  }
  _out.setDevice(&_socket);
  _out.setByteOrder(QDataStream::LittleEndian);
}

void StreamDataWriter::appendPhasedVariantData()
{
  _out << _allele1;
  _out << _allele2;
  if (_printDS) {
    _out << _nAlleles - 1;
    for (int j = 1; j < _nAlleles; ++j) {
      _out << _dose[j];
    }
  }else{
    _out << 0;
  }

  if (_printGP) {
    _out << _gtProbs.length();
    for (int j = 0; j < _gtProbs.length(); ++j) {
      _out << _gtProbs[j];
    }
  }else{
    _out << 0;
  }
}

void StreamDataWriter::finishAndWriteRec()
{
  _out << _marker.chrom();
  _out << _marker.pos();
  _out << _marker.id();
  _out << _nAlleles;
  for(int i=0; i<_nAlleles; i++)
    _out << _marker.allele(i);

  /*
  if (_printDS || _printGP) {
    _out << _r2Est.allelicR2();
    _out << _r2Est.doseR2();
    _out << _nAlleles;
    for (int j = 1; j < _nAlleles; ++j) {
      double af = _cumAlleleProbs[j] / (2 * _r2Est.nGenotypes());
      _out << af;
    }
    _out << _isImputed[_mNum];
  }
  */
  _socket.flush();
}

void StreamDataWriter::finalizeForWindow()
{
  _socket.disconnectFromServer();
}

int main(int argc, char* argv[])
{
  QCoreApplication app(argc, argv);

  // There is no reason why these would fail, but if we can't open up
  // our communication channels, we need to bail and let the parent
  // process consider this a non-starter.

  QFile stderrFile;
  if(!stderrFile.open(stderr, QFile::WriteOnly)){
    return 1;
  }

  ControlStream control(&stderrFile);
  control.sendDebug("startup");

  ImputeOpts opts;
  QString err;
  if (!opts.parseArgs(app.arguments(), err)) {
    control.sendError(err);
    return 1;
  }
  if(opts.pipeName.isEmpty()) {
    control.sendError("No --pipename param provided.\n\nExpected a QLocalServer server name to communicate over with a QDataStream.\n\nIf you are seeing this, you are probably trying to run this program outside of the SVS sub-process context it was designed for.");
    return 1;
  }

  control.sendDebug("args parsed!");

  QThread::msleep(30000);
  StreamTargetDataReader tr(opts.pipeName, control, opts.targetFilePath);

  StreamDataWriter dw(opts.pipeName, control, tr.samples());

  if (opts.readRefData()) {
    // Reference and target data both exist. Open the reference data.
    StreamRefDataReader rr(opts.pipeName, control, opts.refFilePath);

    AllData ad;
    ImputeDriver::phaseAndImpute(ad, tr, rr, dw, opts.window(), opts);

    /*
    if (tr.restrictedZeroMarkerCnt())
      printf(
          "\nWarning: %d positions of the reference window failed to"
          "\n         overlap any target markers.\n",
          tr.restrictedZeroMarkerCnt());

    if (tr.restrictedSingleMarkerCnt())
      printf(
          "\nWarning: %d positions of the reference window overlapped"
          "\n         only one target marker.\n",
          tr.restrictedSingleMarkerCnt());
    */
    return 0;
  } else {
    // Only target data exists. Use the "RefDataReader" base class, which
    // will act as a dummy/placeholder class.
    RefDataReader rr;

    TargetData td;
    ImputeDriver::phaseAndImpute(td, tr, rr, dw, opts.window(), opts);

    return 0;
  }
};

// // For any qt signal/slot MOC stuff
#include "svs_impute_client.moc"

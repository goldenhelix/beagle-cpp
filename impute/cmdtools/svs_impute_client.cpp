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
#include <QTcpSocket>
#include <QThread>
#include <QEventLoop>

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

// Use a QTcpSocket to connect to provided server and read all data
// until close, returning in QByteArray
class SocketReader : public QObject
{
  Q_OBJECT;

public:
  SocketReader(QString pipeName)
    : _pipeName(pipeName)
  {}

  QByteArray readAll(QString &outErr){
    int port = _pipeName.toInt();
    _socket.connectToHost("127.0.0.1", port, QIODevice::ReadOnly);
    if(_socket.state() != QAbstractSocket::ConnectedState){
      QEventLoop wait;
      connect(&_socket, SIGNAL(connected()), &wait, SLOT(quit()));
      wait.exec();
    }
    _buffer.open(QIODevice::WriteOnly);

    connect(&_socket, SIGNAL(disconnected()), this, SLOT(socketDisconnected()));
    connect(&_socket, SIGNAL(readyRead()), this, SLOT(readBuffer()));

    QEventLoop wait;
    connect(this, SIGNAL(readFinished()), &wait, SLOT(quit()));
    wait.exec();

     readBuffer(); // finish read
    _buffer.close();
    return _buffer.data();
  }

signals:
  void readFinished();

public slots:
  void readBuffer(){
    QByteArray read = _socket.readAll();
    //qDebug("read %d", read.size());
    _buffer.write(read);
  }

  void socketDisconnected(){
    //qDebug("disconnected");
    emit readFinished();
  }

private:
  QString _pipeName;
  QTcpSocket _socket;
  QBuffer _buffer;
};

class SocketWriter : public QObject
{
  Q_OBJECT;

public:
  SocketWriter(QTcpSocket* socket)
    : _socket(socket)
  {
    connect(socket, SIGNAL(bytesWritten(qint64)), this, SLOT(markWritten(qint64)));
  }

  void writeBlocking(const QByteArray& data){
    _toWrite = data.size();
    _socket->write(data);
    QEventLoop wait;
    connect(this, SIGNAL(finished()), &wait, SLOT(quit()));
    wait.exec();
  }

signals:
  void finished();

public slots:
  void markWritten(qint64 written){
    _toWrite -= written;
    if(_toWrite <= 0)
      emit finished();
  }
private:
  QTcpSocket* _socket;
  qint64 _toWrite;
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
};

StreamDataParser::StreamDataParser(QString pipeName, ControlStream& control, QString inputToken)
  : _pipeName(pipeName), _control(control), _inputToken(inputToken)
{
}

bool StreamDataParser::readSampleNames()
{
  _control.sendMessage("READ_SAMPLES", _inputToken);

  SocketReader reader(_pipeName);
  QString err;
  QByteArray data = reader.readAll(err);
  if(!err.isEmpty()){
    _control.sendError(err);
    exit(1);
  }

  QDataStream input(data);
  input.setByteOrder(QDataStream::LittleEndian);

  input >> _nSamples;
  if(_nSamples == 0)
    return false;

  for(int i=0; i<_nSamples; i++){
    QByteArray sample;
    input >> sample;
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

  SocketReader reader(_pipeName);
  QString err;
  QByteArray data = reader.readAll(err);
  if(!err.isEmpty()){
    _control.sendError(err);
    exit(1);
  }
  QDataStream input(data);
  input.setByteOrder(QDataStream::LittleEndian);

  int numRecords;
  input >> numRecords;

  // Don't bother asking for more data
  if(numRecords == 0) {
    return false;
  }

  // Should be ready to read numRecords into our buffer
  for(int r=0; r<numRecords; r++){
    Record record;
    input >> record.chrom;
    input >> record.pos;
    input >> record.id;
    input >> record.alleles;
    input >> record.var1;
    input >> record.var2;
    input >> record.phased;

    Q_ASSERT(record.var1.count() > 0 && record.var2.count() > 0 && record.phased.count() > 0);
    Q_ASSERT(record.var1.count() == record.var2.count() && record.var2.count() == record.phased.count());

    if(record.alleles.size() > 0)
      _recordBuffer << record; //Don't add all missing markers
  }
  qDebug("read %d records in %d bytes", numRecords, data.size());
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
  void initializeWindowBuffering(const int initSize, const int nMarkers);
  void appendPhasedVariantData();
  void finishAndWriteRec();
  void finalizeForWindow();

private:
  void outputMarker();
  void outputInfo();
  void outputFormat();

  QString _pipeName;
  ControlStream& _control;
  QBuffer _buffer;
  QDataStream _out;
  QTcpSocket _socket;
};

StreamDataWriter::StreamDataWriter(QString pipeName, ControlStream& control, const Samples& samples)
  : ImputeDataWriter(samples), _pipeName(pipeName), _control(control)
{
}

void StreamDataWriter::writeHeader()
{
  // NO-OP
}

void StreamDataWriter::initializeWindowBuffering(const int initSize, const int nMarkers)
{
  _control.sendDebug(QString("Init window %1 %2").arg(initSize).arg(nMarkers));
  // connect;
  _control.sendMessage("WRITE_RECORDS");

  _buffer.open(QIODevice::WriteOnly | QIODevice::Truncate);
  _out.setDevice(&_buffer);
  _out.setByteOrder(QDataStream::LittleEndian);
  _out << nMarkers;
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
}

void StreamDataWriter::finalizeForWindow()
{
  _buffer.close();
  QByteArray data = _buffer.data();
  qDebug("finished window. Sending %d bytes", data.size());

  _socket.connectToHost("127.0.0.1", _pipeName.toInt(), QIODevice::WriteOnly);
  if(_socket.state() != QAbstractSocket::ConnectedState){
    QEventLoop wait;
    QObject::connect(&_socket, SIGNAL(connected()), &wait, SLOT(quit()));
    wait.exec();
  }
  SocketWriter writer(&_socket);
  writer.writeBlocking(data); //Wait till all bytes written
  _socket.close();
  // qDebug("written!");
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
    control.sendError("No --pipename param provided.\n\nExpected a QTcpServer server name to communicate over with a QDataStream.\n\nIf you are seeing this, you are probably trying to run this program outside of the SVS sub-process context it was designed for.");
    return 1;
  }

  control.sendDebug("args parsed!");

  // QThread::msleep(15000);
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

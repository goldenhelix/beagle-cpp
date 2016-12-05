/* Copyright 2016 Golden Helix, Inc. */
#include <QBuffer>
#include <QCoreApplication>
#include <QDate>
#include <QDir>
#include <QFile>
#include <QStringList>
#include <QTextStream>
#include <QTimer>

#include "impute/cmdtools/imputeopts.h"

// ---------------
// SimpleVcfParser
// ---------------
class SimpleVcfParser
{
public:
  SimpleVcfParser(QString path);
  bool readNextRec();
  QList<QByteArray> sampleNames() const { return _sampleNames; }
  // const char *fileName() { return _file.fileName().toLatin1().constData(); }
  QString fileName();

  // Public data:

  QByteArray curChr;
  int curPos;
  QByteArray curId;
  QList<QByteArray> curAlleles;

  QVector<int> curVar1;
  QVector<int> curVar2;
  QVector<bool> curPhased;

private:
  int _nSamples;
  int _nSamplesP9;
  QList<QByteArray> _sampleNames;
  QFile _file;
};

SimpleVcfParser::SimpleVcfParser(QString path)
{
  _file.setFileName(path);
  if (!_file.open(QFile::ReadOnly)) {
    qWarning("File %s could not be opened!", toC(fileName()));
    exit(1);
  }
  curPos = 0;
}

QString SimpleVcfParser::fileName()
{
  return _file.fileName();
}

bool SimpleVcfParser::readNextRec()
{
  if (!_file.isOpen() || _file.atEnd())
    return false;

  curChr.clear();

  // Clear curChr when it is "consumed", forcing us to scan for next
  while (curChr.isEmpty()) {
    if (_file.atEnd())
      return false;
    QByteArray line = _file.readLine();
    if (line.startsWith("#CHROM")) {
      _sampleNames = line.split('\t').mid(9);

      _nSamples = _sampleNames.length();
      int lastSamp = _nSamples - 1;
      int lastSampLengthM1 = _sampleNames[lastSamp].length() - 1;
      if (_sampleNames[lastSamp][lastSampLengthM1] == '\n')
        _sampleNames[lastSamp] = _sampleNames[lastSamp].mid(0, lastSampLengthM1);

      _nSamplesP9 = _nSamples + 9;

      curVar1.fill(0, _nSamples);  // Allocate to the correct length.
      curVar2.fill(0, _nSamples);
      curPhased.fill(true, _nSamples);
    }
    if (line.startsWith("#"))
      continue;

    QList<QByteArray> fields = line.split('\t');
    if (fields.size() < _nSamplesP9) {
      qWarning("A line in file %s (after record %s/%d) has %d fields instead of %d",
               toC(fileName()), curChr.constData(), curPos, fields.size(), _nSamplesP9);
      exit(1);
    }

    int nSamps = fields.length();
    int lastSamp = nSamps - 1;
    int lastSampLengthM1 = fields[lastSamp].length() - 1;
    if (fields[lastSamp][lastSampLengthM1] == '\n')
      fields[lastSamp] = fields[lastSamp].mid(0, lastSampLengthM1);

    curChr = fields[0];
    curPos = fields[1].toInt();
    curId = fields[2];
    curAlleles.clear();
    curAlleles << fields[3];
    curAlleles << fields[4].split(',');

    curPhased.fill(true);
    for (int fieldNum = 9; fieldNum < _nSamplesP9; fieldNum++) {
      QByteArray thisField = fields[fieldNum].split(':')[0];
      int symbolLoc = thisField.indexOf('|');
      if (symbolLoc < 0) {
        symbolLoc = thisField.indexOf('/');
        if (symbolLoc < 0) {
          qWarning("Bad variant data for record %s/%d in file %s!", curChr.constData(), curPos,
                   toC(fileName()));
          exit(1);
        }

        curPhased[fieldNum - 9] = false;
      }

      QByteArray thisLeft = thisField.mid(0, symbolLoc);
      curVar1[fieldNum - 9] = (thisLeft == ".") ? -1 : thisLeft.toInt();

      QByteArray thisRight = thisField.mid(symbolLoc + 1);
      curVar2[fieldNum - 9] = (thisRight == ".") ? -1 : thisRight.toInt();
    }

    return true;
  }
  return false;
}

// -----------------
// TextRefDataReader
// -----------------
class TextRefDataReader : public RefDataReader
{
public:
  TextRefDataReader(QString path);

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
  SimpleVcfParser _parser;
};

TextRefDataReader::TextRefDataReader(QString path) : _parser(path)
{
  // Cache the first record, assuming it exists. Else mark that we
  // don't have a "next" record. Also set up the "samples" object.

  _nextRecordExists = _parser.readNextRec();

  if (!_nextRecordExists) {
    qWarning("No data exists for file:  %s !", toC(_parser.fileName()));
    exit(1);
  }

  foreach (QByteArray sampName, _parser.sampleNames())
    _samples.setSamp(SampleNames::getIndex(sampName));

  _nSamples = _samples.nSamples();

  assembleNextRefRec();
}

BitSetRefGT TextRefDataReader::nextRec() const
{
  return _nextRecord;
}

void TextRefDataReader::advanceRec()
{
  // Read the record AFTER the "next" record. (The "next" record is
  // the one already in memory.) Mark if it doesn't read, which should
  // usually happen only if we are at the EOF.

  _nextRecordExists = _parser.readNextRec();

  if (_nextRecordExists)
    assembleNextRefRec();
}

void TextRefDataReader::assembleNextRefRec()
{
  BitSetRefGT newrr(_samples);

  newrr.setIdInfo(ChromeIds::getIndex(_parser.curChr), _parser.curPos, _parser.curId);

  foreach (QByteArray al, _parser.curAlleles)
    newrr.addAllele(al);

  // Go through checks to make sure that this really is a reference record....
  for (int samp = 0; samp < _nSamples; ++samp) {
    int a1 = _parser.curVar1[samp];
    int a2 = _parser.curVar2[samp];

    if (!_parser.curPhased[samp]) {
      qWarning("Record %s/%d of file %s is not phased!", _parser.curChr.constData(), _parser.curPos,
               toC(_parser.fileName()));
      exit(1);
    }

    if (a1 == -1) {
      qWarning("Record %s/%d of file %s has missing values!", _parser.curChr.constData(),
               _parser.curPos, toC(_parser.fileName()));
      exit(1);
    }

    if (a2 == -1) {
      qWarning("Record %s/%d of file %s has missing values!", _parser.curChr.constData(),
               _parser.curPos, toC(_parser.fileName()));
      exit(1);
    }
  }

  newrr.storePhasedAlleles(_parser.curVar1, _parser.curVar2);

  _nextRecord = newrr;
}

// -----------------
// TextTargetDataReader
// -----------------
class TextTargetDataReader : public TargDataReader
{
public:
  TextTargetDataReader(QString path);

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
  SimpleVcfParser _parser;
};

TextTargetDataReader::TextTargetDataReader(QString path) : _parser(path)
{
  // Cache the first record, assuming it exists. Else mark that we
  // don't have a "next" record. Also set up the "samples" object.

  _nextRecordExists = _parser.readNextRec();

  if (!_nextRecordExists) {
    qWarning("No data exists for file:  %s !", toC(_parser.fileName()));
    exit(1);
  }

  foreach (QByteArray sampName, _parser.sampleNames())
    _samples.setSamp(SampleNames::getIndex(sampName));

  assembleNextTargetRec();
}

BitSetGT TextTargetDataReader::nextRec() const
{
  return _nextRecord;
}

void TextTargetDataReader::advanceRec()
{
  // Read the record AFTER the "next" record. (The "next" record is
  // the one already in memory.) Mark if it doesn't read, which should
  // usually happen only if we are at the EOF.

  _nextRecordExists = _parser.readNextRec();

  if (_nextRecordExists)
    assembleNextTargetRec();
}

void TextTargetDataReader::assembleNextTargetRec()
{
  BitSetGT newtr(_samples);

  newtr.setIdInfo(ChromeIds::getIndex(_parser.curChr), _parser.curPos, _parser.curId);

  foreach (QByteArray al, _parser.curAlleles)
    newtr.addAllele(al);

  newtr.storeAlleles(_parser.curVar1, _parser.curVar2, _parser.curPhased);

  _nextRecord = newtr;
}

// --------------
// VcfDataWriter
// --------------
class VcfDataWriter : public ImputeDataWriter
{
public:
  VcfDataWriter(QFile& out, const Samples& samples) : ImputeDataWriter(samples), _out(out) {}
  void writeHeader();
  void writeEOF() {}
protected:
  void initializeWindowBuffering(const int initSize);
  void appendPhasedVariantData();
  void finishAndWriteRec();

private:
  void outputMarker();
  void outputInfo();
  void outputFormat();

  QFile& _out;
  QBuffer _recs;
};

static QByteArray asTrimString2(double value)
{
  QByteArray valueString = QByteArray::number(value, 'f', 2);

  int finalPlace = valueString.length() - 1;

  while(finalPlace > 0  &&  valueString[finalPlace] == '0')
  {
    valueString.chop(1);
    finalPlace--;
  }

  if(valueString[finalPlace] == '.')
    valueString.chop(1);

  return valueString;
}

void VcfDataWriter::writeHeader()
{
  _out.write("##fileformat=VCFv4.2\n");
  _out.write("##filedate=" + QDate::currentDate().toString("yyyyMMdd").toLocal8Bit() + "\n");
  _out.write("##source=\"beagle.27Jul16.86a.jar (version 4.1)\"\n");
  _out.write(
      "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated ALT Allele Frequencies\">\n");
  _out.write(
      "##INFO=<ID=AR2,Number=1,Type=Float,Description=\"Allelic R-Squared: estimated squared "
      "correlation between most probable REF dose and true REF dose\">\n");
  _out.write(
      "##INFO=<ID=DR2,Number=1,Type=Float,Description=\"Dosage R-Squared: estimated squared "
      "correlation between estimated REF dose [P(RA) + 2*P(RR)] and true REF dose\">\n");
  _out.write("##INFO=<ID=IMP,Number=0,Type=Flag,Description=\"Imputed marker\">\n");
  _out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
  if (_printDS)
    _out.write(
        "##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"estimated ALT dose [P(RA) + "
        "P(AA)]\">\n");
  if (_printGP)
    _out.write(
        "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Estimated Genotype Probability\">\n");
  _out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

  for (int s = 0, n = _samples.nSamples(); s < n; s++)
    _out.write("\t" + _samples.name(s));
  _out.write("\n");
}

void VcfDataWriter::initializeWindowBuffering(const int initSize)
{
  _recs.close();
  _recs.setData(QByteArray());
  _recs.open(QBuffer::ReadWrite);
}

void VcfDataWriter::appendPhasedVariantData()
{
  _recs.write("\t");
  _recs.write(QByteArray::number(_allele1));
  _recs.write("|");
  _recs.write(QByteArray::number(_allele2));

  if (_printDS) {
    for (int j = 1; j < _nAlleles; ++j) {
      _recs.write((j == 1) ? ":" : ",");
      _recs.write(asTrimString2(_dose[j]));
    }
  }

  if (_printGP) {
    for (int j = 0; j < _gtProbs.length(); ++j) {
      _recs.write((j == 0) ? ":" : ",");
      _recs.write(asTrimString2(_gtProbs[j]));
    }
  }
}

void VcfDataWriter::finishAndWriteRec()
{
  _out.write(_marker.chrom());
  _out.write("\t");
  _out.write(QByteArray::number(_marker.pos()));
  _out.write("\t");
  _out.write(_marker.id());
  _out.write("\t");
  _out.write(_marker.allele(0));
  _out.write("\t");
  _out.write((_marker.nAlleles() == 1) ? "." : _marker.allele(1));
  _out.write("\t.\tPASS\t");  // QUAL  FILTER

  if (_printDS || _printGP) {
    _out.write("AR2=");
    _out.write(QByteArray::number(_r2Est.allelicR2(), 'f', 2));
    _out.write(";DR2=");
    _out.write(QByteArray::number(_r2Est.doseR2(), 'f', 2));

    for (int j = 1; j < _nAlleles; ++j) {
      _out.write((j == 1) ? ";AF=" : ",");
      double af = _cumAlleleProbs[j] / (2 * _r2Est.nGenotypes());
      if(af >= .01  ||  af < .0001)
	_out.write(asTrimString2(af));
      else
	_out.write(QByteArray::number(af, 'f', 4));
    }

    if (_isImputed[_mNum])
      _out.write(";IMP");
  } else {
    _out.write(".");
  }
  _out.write("\t");

  if (_printDS)
    _out.write((_printGP) ? "GT:DS:GP" : "GT:DS");
  else
    _out.write("GT");

  _out.write(_recs.buffer());  // FORMAT
  _out.write("\n");

  _recs.close();
  _recs.setData(QByteArray());
  _recs.open(QBuffer::ReadWrite);
}

int main(int argc, char* argv[])
{
  QCoreApplication app(argc, argv);

  QString err;
  ImputeOpts opts;
  if (!opts.parseArgs(app.arguments(), err)) {
    if(!err.isEmpty()){
      fprintf(stderr, "%s\n\nUse -h for help.\n", toC(err));
    }
    return 1;
  }

  TextTargetDataReader tr(opts.targetFilePath);

  QFile out;
  if (!opts.outFilePath.isEmpty()) {
    out.setFileName(opts.outFilePath);
    if (!out.open(QIODevice::WriteOnly)) {
      printf("Unable to write to %s.\n", opts.outFilePath.toLocal8Bit().constData());
      return 1;
    }
  } else {
    if (!out.open(stdout, QIODevice::WriteOnly)) {
      printf("Unable to open stdout for writing\n");
      return 1;
    }
  }
  VcfDataWriter dw(out, tr.samples());

  if (opts.refFilePath.length()) {
    // Reference and target data both exist. Open the reference data.
    TextRefDataReader rr(opts.refFilePath);

    AllData ad;
    ImputeDriver::phaseAndImpute(ad, tr, rr, dw, opts.window(), opts);

	out.close();

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

    return 0;
  } else {
    // Only target data exists. Use the "RefDataReader" base class, which
    // will act as a dummy/placeholder class.
    RefDataReader rr;

    TargetData td;
    ImputeDriver::phaseAndImpute(td, tr, rr, dw, opts.window(), opts);

    out.close();
    return 0;
  }
};

// // For any qt signal/slot MOC stuff
// #include "imputec.moc"

/* Copyright 2016 Golden Helix, Inc. */
#include <stdio.h>

#include <QCoreApplication>
#include <QDir>
#include <QDate>
#include <QFile>
#include <QStringList>
#include <QTextStream>
#include <QBuffer>
#include <QTimer>

#include "impute/imputedriver.h"

class CmdOpts : public Par
{
public:

  QString outFilePath;
  QString refFilePath;
  QString targetFilePath;

  int window() const {return 4;}
  int overlap() const {return 2;}
  int nThreads() const {return 1;}
  int nSamplingsPerIndividual() const { return 4; }
  int burnin_its() const { return 4; }
  int phase40_its() const { return 4; }
  int niterations() const { return 0; }

  bool parseArgs(QStringList args);

private:
  int _window;
  int _overlap;
  int _nThreads;
  int _nSamplingsPerIndividual;
  int _burnin_its;
  int _phase40_its;
  int _niterations;
};

void printUsage(FILE* fh=stdout)
{
  //fieldList=all is default fieldList=segment means get chr/start/stop only
  QTextStream out(fh);

  out << "Usage: " << qApp->arguments().at(0) << " [params] refpanel.vcf targetpanel.vcf\n\n";
  out << "Params:\n";
  out << "  --out=file_name       Sends output to file (default: stdout)\n";
  out << "  --window=4\n";
  out << "  --overlap=2\n";
  out << "  --nThreads=1\n";
  out << "  --nSamplingsPerIndividual=4\n";
  out << "  --burnin_its=4\n";
  out << "  --phase40_its=4\n";
  out << "  --niterations=0\n";
  out<< "\n";
}

bool CmdOpts::parseArgs(QStringList args)
{
  _window = 4;
  _overlap = 2;
  _nThreads = 1;
  _nSamplingsPerIndividual = 4;
  _burnin_its = 4;
  _phase40_its = 4;
  _niterations = 0;

  args.takeFirst(); //Pop app name

  bool seenSource=false;
  foreach(QString arg, args){
    if(arg == "-h" || arg == "--help"){
      printUsage();
      return false;
    }

    if(arg.startsWith("--")){
      int eqIdx = arg.indexOf("=");
      if(eqIdx < 0) continue;
      bool ok;

      QString opt = arg.mid(2, eqIdx-2);
      QString param = arg.mid(eqIdx+1);
      if(param.startsWith("\"") && param.endsWith("\""))
        param = param.mid(1,param.length()-1);

      if(opt.toLower() == "out"){
        outFilePath = param;
      }
      if(opt.toLower() == "window"){
        _window = param.toInt(&ok);
        if(!ok){
          printf("Param window must be an integer.\n");
          return false;
        }
      }
      if(opt.toLower() == "overlap"){
        _overlap = param.toInt(&ok);
        if(!ok){
          printf("Param overlap must be an integer.\n");
          return false;
        }
      }
      if(opt.toLower() == "nThreads"){
        _nThreads = param.toInt(&ok);
        if(!ok){
          printf("Param nThreads must be an integer.\n");
          return false;
        }
      }
      if(opt.toLower() == "nSamplingsPerIndividual"){
        _nSamplingsPerIndividual = param.toInt(&ok);
        if(!ok){
          printf("Param nSamplingsPerIndividual must be an integer.\n");
          return false;
        }
      }
      if(opt.toLower() == "_burninits"){
        _burnin_its = param.toInt(&ok);
        if(!ok){
          printf("Param burnin_its must be an integer.\n");
          return false;
        }
      }
      if(opt.toLower() == "_phase40its"){
        _phase40_its = param.toInt(&ok);
        if(!ok){
          printf("Param phase40_its must be an integer.\n");
          return false;
        }
      }
      if(opt.toLower() == "niterations"){
        _niterations = param.toInt(&ok);
        if(!ok){
          printf("Param niterations must be an integer.\n");
          return false;
        }
      }
      continue;
    }
    if(refFilePath.isEmpty())
      refFilePath = arg;
    if(targetFilePath.isEmpty())
      targetFilePath = arg;
  }
  if(refFilePath.isEmpty() || targetFilePath.isEmpty()) {
    printUsage(stderr);
    return false;
  }
  return true;
}

// ---------------
// SimpleVcfParser
// ---------------
class SimpleVcfParser {
public:
  SimpleVcfParser(QString path);
  bool hasNextRec() const { return !curChr.isEmpty(); }
  bool readNextRec();
  QList<QByteArray> sampleNames() const { return _sampleNames; }

  QByteArray curChr;
  int curPos;
  QByteArray curId;
  QList<QByteArray> curAlleles;

  QVector<int> curRv1;
  QVector<int> curRv2;

private:
  QList<QByteArray> _sampleNames;
  QFile _file;
};

SimpleVcfParser::SimpleVcfParser(QString path)
{
  _file.setFileName(path);
  _file.open(QFile::ReadOnly);
}

bool SimpleVcfParser::readNextRec()
{
  if(_file.isOpen() || _file.atEnd())
    return false;

  curChr.clear();

  // Clear curChr when it is "consumed", forcing us to scan for next
  while(curChr.isEmpty()){
    if(_file.atEnd())
      return false;
    QByteArray line = _file.readLine();
    if(line.startsWith("#CHROM"))
      _sampleNames = line.split('\t').mid(9);
    if(line.startsWith("#"))
      continue;

    QList<QByteArray> fields = line.split('\t');
    if(fields.size() < 10 + _sampleNames.size()) {
      qWarning("Line with %d fields instead of %d", fields.size(), 10 + _sampleNames.size());
      continue;
    }
    curChr = fields[0];
    curPos = fields[1].toInt();
    curId = fields[2];
    curAlleles.clear();
    curAlleles << fields[3];
    curAlleles << fields[4].split(',');
    // TODO: Parse phased genotypes in fields 10 : 10+_sampleNames.size() and store in curRv1,curRv2
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

  bool canAdvanceWindow() const { return _parser.hasNextRec(); }
  bool hasNextRec() const { return _parser.hasNextRec(); }

  BitSetRefGT nextRec() const;

  void advanceRec();

  // bool lastWindowOnChrom() const;  // Use default implementation here.

private:
  SimpleVcfParser _parser;
};

TextRefDataReader::TextRefDataReader(QString path)
  : _parser(path)
{
  _parser.readNextRec();
}


BitSetRefGT TextRefDataReader::nextRec() const
{
  return BitSetRefGT();
}


void TextRefDataReader::advanceRec()
{
  _parser.readNextRec();
}


// -----------------
// TextRefDataReader
// -----------------
class TextTargetDataReader : public TargDataReader
{
public:
  TextTargetDataReader(QString path);

  bool canAdvanceWindow() const { return _parser.hasNextRec(); }
  bool hasNextRec() const { return _parser.hasNextRec(); }

  BitSetGT nextRec() const;

  void advanceRec();

private:
  SimpleVcfParser _parser;
};

TextTargetDataReader::TextTargetDataReader(QString path)
  : _parser(path)
{
  _parser.readNextRec();
}


BitSetGT TextTargetDataReader::nextRec() const
{
  return BitSetGT();
}

void TextTargetDataReader::advanceRec()
{
  _parser.readNextRec();
}

// --------------
// VcfDataWriter
// --------------
class VcfDataWriter : public ImputeDataWriter
{
public:
  VcfDataWriter(QFile& out, const Samples &samples) : ImputeDataWriter(samples) , _out(out) {}

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

void VcfDataWriter::writeHeader()
{
  _out.write("##fileformat=VCFv4.2\n");
  _out.write("##filedate=" + QDate::currentDate().toString("yyyyMMdd").toLocal8Bit() + "\n");
  _out.write("##source=\"beagle.27Jul16.86a.jar (version 4.1)\"\n");
  _out.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated ALT Allele Frequencies\">\n");
  _out.write("##INFO=<ID=AR2,Number=1,Type=Float,Description=\"Allelic R-Squared: estimated squared correlation between most probable REF dose and true REF dose\">\n");
  _out.write("##INFO=<ID=DR2,Number=1,Type=Float,Description=\"Dosage R-Squared: estimated squared correlation between estimated REF dose [P(RA) + 2*P(RR)] and true REF dose\">\n");
  _out.write("##INFO=<ID=IMP,Number=0,Type=Flag,Description=\"Imputed marker\">\n");
  _out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
  if (_printDS)
    _out.write("##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"estimated ALT dose [P(RA) + P(AA)]\">\n");
  if(_printGP)
    _out.write("##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Estimated Genotype Probability\">\n");
  _out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

  for(int s=0, n=_samples.nSamples(); s<n; s++)
    _out.write("\t" + _samples.name(s));
  _out.write("\n");
}

void VcfDataWriter::initializeWindowBuffering(const int initSize)
{
  _recs.close();
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
      _recs.write(QByteArray::number(_dose[j], 'f', 2));
    }
  }

  if (_printGP) {
    for (int j = 0; j < _gtProbs.length(); ++j) {
      _recs.write((j == 0) ? ":" : ",");
      _recs.write(QByteArray::number(_gtProbs[j], 'f', 2));
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
  _out.write((_marker.nAlleles() == 1) ? "?" : _marker.allele(1) );
  _out.write("\t?\tPASS\t");// QUAL  FILTER

  if (_printDS || _printGP) {
    _out.write("AR2=");
    _out.write(QByteArray::number(_r2Est.allelicR2(), 'f', 2));
    _out.write(";DR2=");
    _out.write(QByteArray::number(_r2Est.doseR2(), 'f', 2));

    for (int j = 1; j < _nAlleles; ++j) {
      _out.write((j == 1) ? ";AF=" : ",");
      _out.write(QByteArray::number(_cumAlleleProbs[j] / (2 * _r2Est.nGenotypes()), 'f', 2));
    }

    if (_isImputed[_mNum])
      _out.write(";IMP");
  } else {
    _out.write("?");
  }
  _out.write("\t");

  if (_printDS)
    _out.write((_printGP) ? "GT:DS:GP" : "GT:DS");
  else
    _out.write("GT");
  _out.write(_recs.buffer());  // FORMAT

  _recs.close();
  _recs.open(QBuffer::ReadWrite);  
}

int main(int argc, char* argv[])
{
  QCoreApplication app(argc,argv);

  CmdOpts opts;
  if(!opts.parseArgs(app.arguments()))
    return 1;


  TextRefDataReader rr(opts.refFilePath);
  TextTargetDataReader tr(opts.targetFilePath);

  QFile out;
  if(!opts.outFilePath.isEmpty()) {
    out.setFileName(opts.outFilePath);
    if(!out.open(QIODevice::WriteOnly)){
      printf("Unable to write to %s.\n", opts.outFilePath.toLocal8Bit().constData());
      return 1;
    }
  }else{
    if(!out.open(stdout, QIODevice::WriteOnly)){
      printf("Unable to open stdout for writing\n");
      return 1;
    }
  }
  VcfDataWriter dw(out, tr.samples());

  TargetData td;
  ImputeDriver::phaseAndImpute(td, tr, rr, dw, opts.window(), opts);

  out.close();
  return 0;
};

// For any qt signal/slot MOC stuff
#include "imputec.moc"

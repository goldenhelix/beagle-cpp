#ifndef IMPUTEOPTS_H
#define IMPUTEOPTS_H
#include <stdio.h>

#include "impute/imputedriver.h"

class ImputeOpts : public Par
{
public:
  QString outFilePath;
  QString refFilePath;
  QString targetFilePath;

  QString pipeName; //Only used in an IPC context

  int window() const { return _window; }    // Test value could be 4 .
  int overlap() const { return _overlap; }  // Test value could be 2 .
  int nThreads() const { return _nThreads; }    // Test value could be 4 .
  int nSamplingsPerIndividual() const
  {
    return _nSamplingsPerIndividual;
  }                                                 // Test value could be 4 .
  bool lowMem() const { return _lowmem; }
  int burnin_its() const { return _burnin_its; }    // Test value could be 4 .
  int phase40_its() const { return _phase40_its; }  // Test value could be 4 .
  int niterations() const { return _niterations; }  // Test value could be 0 .
  bool impute() const { return _impute; }           // Test value could be true .
  bool gprobs() const { return _gprobs; }
  float cluster() const { return _cluster; }
  float ne() const { return _ne; }
  float err() const { return _err; }
  bool readRefData() const { return !refFilePath.isEmpty(); }
  bool parseArgs(QStringList args, QString& outErr);

private:
  int _window;
  int _overlap;
  int _nThreads;
  int _nSamplingsPerIndividual;
  bool _lowmem;
  int _burnin_its;
  int _phase40_its;
  int _niterations;
  bool _impute;
  float _cluster;
  float _ne;
  float _err;
  bool _gprobs;
};

void printUsage(FILE* fh = stdout)
{
  // fieldList=all is default fieldList=segment means get chr/start/stop only
  QTextStream out(fh);

  out << "\n";
  out << "Usage: " << qApp->arguments().at(0) << " [params] refpanel.vcf targetpanel.vcf\n";
  out << "         or\n";
  out << "       " << qApp->arguments().at(0) << " [params] targetpanel.vcf\n\n";
  out << "Params:\n";
  out << "  --out=file_name       Sends output to file (default: stdout)\n";
  out << "  --window=50000\n";
  out << "  --overlap=3000\n";
  out << "  --nThreads=4\n";
  out << "  --nSamplingsPerIndividual=4\n";
  out << "  --lowmem=true\n";
  out << "  --burninits=5\n";
  out << "  --phase40its=5\n";
  out << "  --niterations=5\n";
  out << "  --impute=true\n";
  out << "  --cluster=0.005\n";
  out << "  --ne=1000000.0\n";
  out << "  --err=0.0001\n";
  out << "  --gprobs=false\n";
  out << "\n";
}

bool ImputeOpts::parseArgs(QStringList args, QString& outErr)
{
  _window = Par::window();
  _overlap = Par::overlap();
  _nThreads = Par::nThreads();
  _nSamplingsPerIndividual = Par::nSamplingsPerIndividual();
  _lowmem = Par::lowMem();
  _burnin_its = Par::burnin_its();
  _phase40_its = Par::phase40_its();
  _niterations = Par::niterations();
  _impute = Par::impute();
  _cluster = Par::cluster();
  _ne = Par::ne();
  _err = Par::err();
  _gprobs = false;

  args.takeFirst();  // Pop app name

  bool seenSource = false;
  foreach (QString arg, args) {
    if (arg == "-h" || arg == "--help") {
      printUsage();
      return false;
    }

    if (arg.startsWith("--")) {
      int eqIdx = arg.indexOf("=");
      if (eqIdx < 0)
        continue;
      bool ok;

      QString opt = arg.mid(2, eqIdx - 2);
      QString param = arg.mid(eqIdx + 1);
      if (param.startsWith("\"") && param.endsWith("\""))
        param = param.mid(1, param.length() - 1);

      if (opt.toLower() == "out") {
        outFilePath = param;
      } else if (opt.toLower() == "pipename") {
        pipeName = param;
      } else if (opt.toLower() == "window") {
        _window = param.toInt(&ok);
        if (!ok) {
          outErr = "Param window must be an integer.";
          return false;
        }
      } else if (opt.toLower() == "overlap") {
        _overlap = param.toInt(&ok);
        if (!ok) {
          outErr = "Param overlap must be an integer.";
          return false;
        }
      }
      else if (opt.toLower() == "nthreads") {
        _nThreads = param.toInt(&ok);
        if (!ok) {
          outErr = "Param nthreads must be an integer.";
          return false;
        }
      }
      else if (opt.toLower() == "nsamplingsperindividual") {
        _nSamplingsPerIndividual = param.toInt(&ok);
        if (!ok) {
          outErr = "Param nsamplingsperindividual must be an integer.";
          return false;
        }
      } else if (opt.toLower() == "lowmem") {
        if (param.toLower() == "true")
          _lowmem = true;
        else if (param.toLower() == "false") {
          _lowmem = false;
        } else {
          outErr = "Param lowmem must be either \"true\" or \"false\".";
          return false;
        }
      } else if (opt.toLower() == "burninits") {
        _burnin_its = param.toInt(&ok);
        if (!ok) {
          outErr = "Param burninits must be an integer.";
          return false;
        }
      } else if (opt.toLower() == "phase40its") {
        _phase40_its = param.toInt(&ok);
        if (!ok) {
          outErr = "Param phase40its must be an integer.";
          return false;
        }
      } else if (opt.toLower() == "niterations") {
        _niterations = param.toInt(&ok);
        if (!ok) {
          outErr = "Param niterations must be an integer.";
          return false;
        }
      } else if (opt.toLower() == "impute") {
        if (param.toLower() == "true")
          _impute = true;
        else if (param.toLower() == "false") {
          _impute = false;
        } else {
          outErr = "Param impute must be either \"true\" or \"false\".";
          return false;
        }
      } else if (opt.toLower() == "cluster") {
        _cluster = param.toFloat(&ok);
        if (!ok) {
          outErr = "Param cluster must be a float.";
          return false;
        }
      } else if (opt.toLower() == "ne") {
        _ne = param.toFloat(&ok);
        if (!ok) {
          outErr = "Param ne must be a float.";
          return false;
        }
      } else if (opt.toLower() == "err") {
        _err = param.toFloat(&ok);
        if (!ok) {
          outErr = "Param err must be a float.";
          return false;
        }
      } else if (opt.toLower() == "gprobs") {
        if (param.toLower() == "true")
          _gprobs = true;
        else if (param.toLower() == "false") {
          _gprobs = false;
        } else {
          outErr = "Param lowmem must be either \"true\" or \"false\".";
          return false;
        }
      } else {
        outErr = "Bad parameter %s.";
        return false;
      }
      continue;
    }
    if (refFilePath.isEmpty())
      refFilePath = arg;
    else if (targetFilePath.isEmpty())
      targetFilePath = arg;
  }

  if (refFilePath.isEmpty() && targetFilePath.isEmpty()) {
    outErr = "No refFilePath or targetFilePath specified.";
    return false;
  }

  if (targetFilePath.isEmpty()) {
    // Target file only.
    targetFilePath = refFilePath;
    refFilePath.clear();
  }
  return true;
}

#endif // IMPUTEOPTS_H

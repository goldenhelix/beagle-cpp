/* Copyright 2015 Golden Helix, Inc. */
#ifndef CSTRING_H
#define CSTRING_H

#include <QByteArray>
#include <QPair>
#include <QString>

// TODO: This class is going to break inheritance from QByteArray and
// grow ability to easily intern char* data allocated outside it's
// scope. It will also grow the ability to refcount a shared char*
// allocation while having it's individualized offset pointer into a
// null-separated list of strings.

//Store UTF8 string data in data arrays
class CString : public QByteArray
{
public:
  CString() : QByteArray() {}
  CString( const char * str )  : QByteArray( str ) {}
  CString( const char * data, int size ) : QByteArray(data, size) {}
  CString( const QByteArray & other ) : QByteArray(other) {}
  CString( const CString & other ) : QByteArray(other) {}
  CString( const QString & other ) : QByteArray(other.toUtf8()) {}
  CString(int size, char ch ) : QByteArray(size, ch) {}

  QString asQString() const { return QString::fromUtf8( constData(), size() ); }
  operator QString () const { return asQString(); }

  //replaces contents with null terminated input string
  void read(const char * str);

  //conversion functions that return missing values (missing.h) when
  //failed to convert


  static const QVector< QPair<CString, CString> >& trueFalsePairVector();
  //true if value could be converted to a boolean
  bool isBool(const QVector< QPair<CString, CString> >* trueFalsePairs = 0) const;
  //one of known "true" strings, otherwise false
  bool toPureBool(const QVector< QPair<CString, CString> >* trueFalsePairs = 0) const;
  //one of known "true" or "false" strings, otherwise missing
  char toBool(const QVector< QPair<CString, CString> >* trueFalsePairs = 0) const;

  int toInt(bool* ok=0) const;
  qint64 toInt64(bool* ok=0) const;
  float toFloat(bool* ok=0) const;
  double toDouble(bool* ok=0) const;

  //Note that you just call constData() to get a const char *, but in
  //most cases it will auto-cast for you.
};

/// Helper macro for sending a QString to c code that requires const char * (equivelant to qPrintable)
#define toC(string) (string).toLocal8Bit().constData()
/// Helper macro for interfacing with any C library. We always use UTF8 on persistant data (and SQLite)
#define to8(string) (string).toUtf8().constData()

#endif

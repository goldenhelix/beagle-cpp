#include "impute/genotypecorrection.h"
#include "impute/haplotypepair.h"

class Edit
{
public:

  /**
   * Constructs a new {@code Edit} instance.
   *
   * @param hapPair a haplotype pair
   * @param marker the marker index
   * @param newAllele1 the post-edit first allele
   * @param newAllele2 the post-edit second allele
   */
  Edit(const HapPair &hapPair, int marker, int newAllele1, int newAllele2);
  Edit(const Edit &other);

  const HapPair &hapPair() const { return _hapPair; }
  int marker() const { return _marker; }
  int newAllele1() const { return _newAllele1; }
  int newAllele2() const { return _newAllele2; }

private:
  const HapPair &_hapPair;
  int _marker;
  int _newAllele1;
  int _newAllele2;
};

Edit::Edit(const HapPair &hapPair, int marker, int newAllele1, int newAllele2)
  : _hapPair(hapPair), _marker(marker), _newAllele1(newAllele1), _newAllele2(newAllele2)
{
  Q_ASSERT_X(marker >= 0  &&  marker < hapPair.nMarkers(),
             "Edit::Edit",
             "marker<0 || marker >= hapPair.nMarkers()");

  Q_ASSERT_X(newAllele1 >= 0  &&  newAllele1 < hapPair.marker(marker).nAlleles(),
             "Edit::Edit",
             "newAllele1<0 || newAllele1 >= hapPair.marker(marker).nAlleles()");

  Q_ASSERT_X(newAllele2 >= 0  &&  newAllele2 < hapPair.marker(marker).nAlleles(),
             "Edit::Edit",
             "newAllele2<0 || newAllele2 >= hapPair.marker(marker).nAlleles()");
}

Edit::Edit(const Edit &other)
  : _hapPair(other._hapPair), _marker(other._marker), _newAllele1(other._newAllele1),
    _newAllele2(other._newAllele2)
{ }

static void copyAlleles(const HapPair &hapPair, QList<int> &alleles1,
                        QList<int> &alleles2)
{
  alleles1.clear();
  alleles2.clear();

  for (int m=0, n=hapPair.nMarkers(); m<n; ++m)
  {
    alleles1.append(hapPair.allele1(m));
    alleles2.append(hapPair.allele2(m));
  }
}

static void checkMarkersAndSamples(const HapPair &hapPair, const MaskedEndsGL &gl)
{
  Q_ASSERT_X(hapPair.markers() == gl.markers(),
             "checkMarkersAndSamples (genotypecorrection.cpp)",
             "hapPair.markers() != gl.markers()");

  Q_ASSERT_X(hapPair.samples() == gl.samples(),
             "checkMarkersAndSamples (genotypecorrection.cpp)",
             "hapPair.samples() != gl.samples()");
}

static QList<Edit> getEdits(const HapPair &hapPair, const MaskedEndsGL &gl /*, Random random */ )
{
  QList<Edit> corrections;
  int sample = hapPair.sampleIndex();
  for (int marker=0, n=gl.nMarkers(); marker<n; ++marker)
  {
    int hapPairA1 = hapPair.allele1(marker);
    int hapPairA2 = hapPair.allele2(marker);
    if (gl.gl(marker, sample, hapPairA1, hapPairA2) <= 0.0f)
    {
      int glA1 = gl.allele1(marker, sample);
      int glA2 = gl.allele2(marker, sample);
      if (glA1 >= 0  &&  glA2 >= 0)
      {

#ifdef SIMULATE_RANDOM
        if (gl.isPhased(marker, sample) == false  &&  (marker % 2)) {
          int tmp = glA1;
          glA1 = glA2;
          glA2 = tmp;
        }
#else
        if (gl.isPhased(marker, sample) == false  &&  (qrand() % 2)) {
          int tmp = glA1;
          glA1 = glA2;
          glA2 = tmp;
        }
#endif

        corrections.append(Edit(hapPair, marker, glA1, glA2));
      }
    }
  }
  return corrections;
}

static HapPair updatedHapPair(const HapPair &hapPair,
                              const QList<Edit> &edits)
{
  if (edits.isEmpty())
    return hapPair;
  else
  {
    QList<int> alleles1;
    QList<int> alleles2;
    copyAlleles(hapPair, alleles1, alleles2);

    for (int j=0, n=edits.size(); j<n; ++j)
    {
      int m = edits[j].marker();
      alleles1[m] = edits[j].newAllele1();
      alleles2[m] = edits[j].newAllele2();
    }

    return HapPair(hapPair.markers(), hapPair.samples(),
                   hapPair.sampleIndex(), alleles1, alleles2);
  }
}

/*
    private static final String headerLine = "MARKER" + Const.tab
                        + "SAMPLE" + Const.tab
                        + "REF" + Const.tab + "ALT" + Const.tab
                        + "INPUT_GT" + Const.tab + "ESTIMATED_GT";
     (After the call to getEdits in run/correct:)
     for (Edit edit : edits) {
         out.println(edit);
     }
*/

void GenotypeCorrection::correct(QList<HapPair> &hapPairs, const MaskedEndsGL &gl, int seed)
{
  qsrand(seed);

  for (int j=0, n=hapPairs.size(); j<n; ++j)
  {
    HapPair hapPair = hapPairs[j];
    checkMarkersAndSamples(hapPair, gl);
    QList<Edit> edits = getEdits(hapPair, gl /*, random */ );
    HapPair revHapPair = updatedHapPair(hapPair, edits);
    hapPairs[j] = revHapPair;
  }
}


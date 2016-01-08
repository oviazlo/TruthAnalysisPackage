/// WprimeHist.h

/// David Adams
/// May 2011

/// Class to create a standard histogram for the Wprime analysis.

#include <string>
#include <iostream>

class TH1;

class WprimeHist {

public:

  typedef std::string Name;

public:

  /// Create from name, title and binning.
  /// These just call the similar TH1F constructors.
  static TH1* create(Name name, Name title, int xbin, double xmin, 
                     double xmax);
  static TH1* create(Name name, Name title, int xbin, double* bins);

  /// Create histogram with log binning.
  static TH1* createlog(Name name, Name title, int xbin, double xmin, 
                        double xmax);

  /// Create a standard Wprime histogram with a name and title.
  /// Return null if the name is not recogized.
  ///   name =  mt, pt, met, eta, ....
  ///   hprefix = prefix for title, e.g. "My sample" --> "My sample pT"
  ///       sel = selection level:
  ///               ""     - default, typically weighted after final selection
  ///               finalu - unweighted after final selection
  ///               genu   - unweighted afte generator selection
  ///               genw   - weighted afte generator selection
  ///             Other values are allowed.
  ///   hname = name assigned to histogram, e.g. "myhist" --> "hmyhist_pt"
  static TH1* standard(Name name, Name hpname, Name sel, Name hptitle);

  /// Create a standard histogram for a sample.
  ///        name =  mt, pt, met, eta, ....
  ///  samplename = wenu, ..., see WprimeSample
  ///         sel = selection level (see above)
  static TH1* sample(Name name, Name samplename, Name sel ="");
  
  /// Get/Set the lepton to use for plot labels.
  /// Returns the old value.
  ///   blank for no change (i.e. when fetching value)
  ///   e for "e"
  ///   mu or #mu for "#mu"
  ///   anything else for "l"
  static Name lepton(Name slepton ="");
};

// WprimeSample.h

// David Adams
// May 2011

// Class to describe each of the wprime samples.

#include <string>
#include <iostream>

class WprimeSample {

public:

  typedef std::string Name;

public:

  // Construct a sample from an alias.
  // If the alias is know, the sample properties (including
  // the name) will be set.
  explicit WprimeSample(Name alias);

  // Return the alias used to construct the sample.
  Name alias() const { return m_alias; }

  // Return the name used to identify the sample.
  Name name() const { return m_name; }

  // Return a description of the sample.
  Name title() const { return m_title; }

  // Return a latex name for the sample.
  Name latex() const { return m_latex; }

  // Return the generator cross section in pb
  // (including BR and generator efficiency).
  double genxspb() const { return m_genxspb; }

  // Return the calculated, typically NNLO or NLO,
  // cross section.
  double calxspb() const { return m_calxspb; }

  // Return the run number.
  int run() const { return m_run; }

private:

  Name m_alias;
  Name m_name;
  Name m_title;
  Name m_latex;
  double m_genxspb;
  double m_calxspb;
  int m_run;

};

// Output stream.
std::ostream& operator<<(std::ostream& lhs, const WprimeSample& rhs);


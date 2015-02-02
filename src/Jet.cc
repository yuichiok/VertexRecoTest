
#include "Jet.hh"
using EVENT::ReconstructedParticle;
using std::vector;
namespace TTbarAnalysis
{
	Jet::Jet (float btag,float ctag, int number,const double * momentum)
	{
		myBTag = btag;
		myCTag = ctag;
		myNumber = number;
		myMomentum = momentum;
		myMCPDG = 0;
		myParticles = NULL;
	}
	float Jet::GetBTag()
	{
		return myBTag;
	}
	float Jet::GetCTag()
	{
		return myCTag;
	}
	int Jet::GetNumberOfVertices()
	{
		return myNumber;
	}
	const double * Jet::GetMomentum()
	{
		return myMomentum;
	}
	const vector< ReconstructedParticle * > * Jet::GetParticles() const
	{
		return myParticles;
	}
	void Jet::SetParticles(const vector< ReconstructedParticle * > & particles)
	{
		myParticles = new vector< ReconstructedParticle * >( particles);
	}
	
	int Jet::GetMCPDG()
	{
		return myMCPDG;
	}
	void Jet::SetMCPDG(int pdg)
	{
		myMCPDG = pdg;
	}

}

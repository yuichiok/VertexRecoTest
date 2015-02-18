#include "Particle.hh"
using std::vector;
using std::string;
using EVENT::Vertex;
using EVENT::ReconstructedParticle;
namespace TTbarAnalysis
{
	Particle::Particle ()
	{
		myMomentum = NULL;
		myVertex = NULL;
	}
	double * Particle::GetMomentum()
	{
		return myMomentum;
	}
	void Particle::SetMomentum(const double * momentum)
	{
		myMomentum = new double[3];
		for (int i = 0; i < 3; i++) 
		{
			myMomentum[i] = momentum[i];
		}
	}
	Vertex * Particle::GetVertex()
	{
		return myVertex;
	}
	float Particle::GetOffset()
	{
		return myOffset;
	}
	void Particle::SetOffset(float offset)
	{
		myOffset = offset;
	}
	void Particle::SetVertex(Vertex * vertex)
	{
		myVertex = vertex;
	}
	float Particle::GetTheta()
	{
		return myTheta;
	}
	void Particle::SetTheta(float theta)
	{
		myTheta = theta;
	}
	
	void Particle::Assign(ReconstructedParticle * initial)
	{
		myMomentum = new double( *initial->getMomentum());
	}
}

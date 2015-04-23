#include "ParticleOperator.hh"
using std::vector;
using std::string;
using EVENT::LCCollection;
using EVENT::MCParticle;
using EVENT::Vertex;
using EVENT::ReconstructedParticle;
namespace TTbarAnalysis 
{
	ParticleOperator:: ParticleOperator()
	{
	}
	bool ParticleOperator::CompareParticles(ReconstructedParticle * particle1, ReconstructedParticle * particle2)
	{
		if (particle1 == particle2) 
		{
			std::cout << "Equal by pointer\n";
			return true;
		}
		//std::cout << "Comparing particles with modules " << recomodule << " & " << mcmodule << "\n";
		if (particle1->getCharge() * particle2->getCharge() < 0.0) 
		{
			//std::cout << "Rejected by charge\n";
			return false;
		}
		float angle = MathOperator::getAngleBtw(particle1->getMomentum(), particle2->getMomentum());
		if (angle > 0.02) 
		{
			//std::cout << "Rejected by angle " << angle << "\n";
			return false;
		}
		float recomodule = MathOperator::getModule(particle2->getMomentum());
		float mcmodule = MathOperator::getModule(particle1->getMomentum());
		float ratio = (1 - recomodule/mcmodule > 0.0) ? 1 - recomodule/mcmodule : recomodule/mcmodule - 1.0;
		if (ratio > 0.04) 
		{
			//std::cout << "Rejected by ratio " << ratio << "\n";
			return false;
		}
		return true;
	}
	bool ParticleOperator::CheckForVertex(EVENT::ReconstructedParticle * particle1)
	{
		Vertex * vertex = particle1->getStartVertex ();
		if (vertex && !vertex->isPrimary()) 
		{
			return true;
		}
		return false;
	}
	bool ParticleOperator::IsDublicate(MCParticle * particle, vector< MCParticle * > & list)
	{
		for (unsigned int i = 0; i < list.size(); i++) 
		{
			if (particle == list[i]) 
			{
				return true;
			}
		}
		return false;
	}
}

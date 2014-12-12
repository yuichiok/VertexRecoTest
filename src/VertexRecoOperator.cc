#include "VertexRecoOperator.hh"
using std::vector;
using std::string;
using EVENT::LCCollection;
using EVENT::Vertex;
using EVENT::ReconstructedParticle;
namespace TTbarAnalysis 
{
	VertexRecoOperator:: VertexRecoOperator() 
	{
		myAngleCut = 0.08;
		myPrecisionCut = 1.5;
	}
	vector< VertexTag * > * VertexRecoOperator::Compare(LCCollection * reconstructed, LCCollection * mc)
	{
		int number = reconstructed->getNumberOfElements();
		vector< VertexTag * > * result = new vector< VertexTag * >();
		for (int i = 0; i < number; i++) 
		{
			Vertex * vertex = dynamic_cast< Vertex * >( reconstructed->getElementAt(i) ) ;
			VertexTag * tagged = Tag(vertex, mc);
			if (tagged) 
			{
				result->push_back(tagged);
			}
			else 
			{
				myUnknownVertexes.push_back(vertex);
			}
		}

		return result;
	}

	vector< VertexTag * > * VertexRecoOperator::CompareDirection(LCCollection * reconstructed, LCCollection * mc)
	{
		int number = reconstructed->getNumberOfElements();
		int mcnumber = mc->getNumberOfElements();
		vector< VertexTag * > * result = new vector< VertexTag * >();
		for (int i = 0; i < number; i++) 
		{
			Vertex * recovertex = dynamic_cast< Vertex * >( reconstructed->getElementAt(i) ) ;
			//vector<float> angles = getAngles(recovertex);
			bool passed = false;
			for (int j = 0; j < mcnumber; j++) 
			{
				Vertex * mcvertex = dynamic_cast< Vertex * >( mc->getElementAt(j) ) ;
				passed = false;
				//vector<float> mcangles = getAngles(mcvertex);
				/*for (int k = 0; k < mcangles.size(); k++) 
				{
					std::cout << "Angle difference: " << mcangles[k] - angles[k] << '\n';
					if (abs(mcangles[k] - angles[k]) > myAngleCut) 
					{
						failed = true;
						//break;
					}
				}*/
				float angle = MathOperator::getAngle(mcvertex->getAssociatedParticle()->getMomentum(), recovertex->getAssociatedParticle()->getMomentum());
				std::cout << "Angle: " << angle  << '\n';
				if (angle < myAngleCut) 
				{
					std::cout << "Vertex tagged with pdg " << mcvertex->getParameters()[1] << '\n';
					passed = true;
					result->push_back(new VertexTag(recovertex, mcvertex));
					break;
				}
			}
			if (!passed) 
			{
				myUnknownVertexes.push_back(recovertex);
			}
		}
		return result;
	}

	VertexTag * VertexRecoOperator::Tag(Vertex * vertex, LCCollection * mc)
	{
		int number = mc->getNumberOfElements();
		VertexTag * result = NULL ;
		float min = 1000.0;
		int candidateNumber = -1; 
		for (int i = 0; i < number; i++) 
		{
			Vertex * mcvertex = dynamic_cast< Vertex * >( mc->getElementAt(i) ) ;
			float distance = MathOperator::getDistance(vertex->getPosition(), mcvertex->getPosition());
			myDistances.push_back(distance);
			if (distance < min && !CheckTaken(mcvertex)) 
			{
				min = distance;
				candidateNumber = i;
			}
		}
		if (min < myPrecisionCut) 
		{
			Vertex * mcvertex = dynamic_cast< Vertex * >( mc->getElementAt(candidateNumber) ) ;
			myTakenVertices.push_back(GetID(mcvertex));

			std::cout << "Found "<< mcvertex->getParameters()[2] << " vertex in " << min << "mm radius and mcid " << GetID(mcvertex) << "\n";
			result = new VertexTag(vertex, mcvertex);
			SetMinimalDistance(result, mc);
		}
		return result;
	}
	void VertexRecoOperator::SetMinimalDistance(VertexTag * tag, LCCollection * mc)
	{
		int number = mc->getNumberOfElements();
		float min = 1000.0;
		Vertex * mcvertex = tag->GetVertex();
		for (int j = 0; j < number; j++) 
		{
			Vertex * mcvertex2 = dynamic_cast< Vertex * >( mc->getElementAt(j) ) ;
			if (((int)mcvertex2->getParameters()[1]) != tag->GetInitialPDG() || ((int)mcvertex2->getParameters()[2]) != tag->GetGeneration()) 
			{
				float distance = MathOperator::getDistance(mcvertex->getPosition(), mcvertex2->getPosition());
				if (distance < min) 
				{
					min = distance;
				}
			}
		}
		std::cout << "Minimal distance "<< min << "\n";
		tag->SetMinimalDistance(min);
	}
	bool VertexRecoOperator::CheckTaken(Vertex * mcVertex)
	{
		bool result = false;
		for (int i = 0; i < myTakenVertices.size(); i++) 
		{
			int id = GetID(mcVertex);
			if (myTakenVertices[i] == id) 
			{
				result = true;
				std::cout << "Vertex is taken!\n";
				break;
			}
		}
		return result;
	}
	vector< float >  VertexRecoOperator::GetDistances()
	{
		return myDistances;
	}
	vector< Vertex * > & VertexRecoOperator::GetUnknownVertexes()
	{
		return myUnknownVertexes;
	}
	int VertexRecoOperator::GetID(EVENT::Vertex * mcVertex)
	{
		return ((int)mcVertex->getParameters()[2])*((int)mcVertex->getParameters()[1]);
	}
	vector<float> VertexRecoOperator::getAngles(EVENT::Vertex * vertex)
	{
		ReconstructedParticle * particle = vertex->getAssociatedParticle();
		const float * position = vertex->getPosition();
		double convert[3];
		for (int i = 0; i < 3; i++) 
		{
			convert[i] = position[i];
		}
		vector<float> direction = MathOperator::getDirection(convert);//particle->getMomentum());
		return MathOperator::getAngles(direction);
	}
} /* TTbarAnalysis */

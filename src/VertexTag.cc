#include "VertexTag.hh"
using std::vector;
using std::string;
using EVENT::Vertex;
namespace TTbarAnalysis 
{
	VertexTag:: VertexTag(Vertex * recoVertex, Vertex * mcVertex)
	{
		myRecoVertex = recoVertex;
		myMCVertex = mcVertex;
		myInitialPDG = myMCVertex->getParameters()[1];
	}
	Vertex * VertexTag::GetVertex()
	{
		return myRecoVertex;
	}
	Vertex * VertexTag::__GetMCVertex()
	{
		return myMCVertex;
	}
	int VertexTag::GetGeneration()
	{
		return myMCVertex->getParameters()[2];
	}
	int VertexTag::GetInitialPDG()
	{
		return myInitialPDG;
	}
	float VertexTag::GetMinimalDistance()
	{
		return myMinimalDistanceMC;
	}
	void VertexTag::SetMinimalDistance(float d)
	{
		myMinimalDistanceMC = d;
	}
	int VertexTag::__GetMCTrackNumber()
	{
		return myMCVertex->getAssociatedParticle()->getParticles().size();
	}
	int VertexTag::GetTrackNumber()
	{
		return myRecoVertex->getAssociatedParticle()->getParticles().size();
	}


} /* TTbarAnalysis */

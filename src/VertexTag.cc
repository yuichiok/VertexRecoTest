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


} /* TTbarAnalysis */

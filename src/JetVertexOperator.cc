#include "JetVertexOperator.hh"
using std::vector;
using std::string;
using EVENT::LCCollection;
using EVENT::MCParticle;
using EVENT::Vertex;
using EVENT::ReconstructedParticle;
namespace TTbarAnalysis 
{
	JetVertexOperator:: JetVertexOperator (LCCollection * pfo, LCCollection * egprongs, LCCollection * trackrel)
	{
		myPFO = pfo;
		myEGProngs = egprongs;
		if (trackrel) 
		{
			myTrackRel = trackrel;
		}
		myAlgorithmName = "lcfiplus";
		myAngleCut = 0.4; // CHANGED
		ip[0] = 0.0;
		ip[1] = 0.0;
		ip[2] = 0.0;

	}
	vector< ReconstructedParticle * > * JetVertexOperator::GetRecoParticles(LCCollection * prongs, LCCollection * rel)
	{
		vector<ReconstructedParticle*> * result = new vector<ReconstructedParticle*>();
		int prongnumber = prongs->getNumberOfElements();
		LCRelationNavigator navigator(rel);
		for (int i = 0; i < prongnumber; i++) 
		{
			MCParticle * particle =  dynamic_cast< MCParticle * >(prongs->getElementAt(i));
			vector< LCObject * > obj = navigator.getRelatedFromObjects(particle);
			if (obj.size() < 1) 
			{
				continue;
			}
			if (obj.size() > 0) 
			{
				float maxweight = 0.0;
				vector< float > weights = navigator.getRelatedFromWeights (particle); 
				ReconstructedParticle * winner = NULL; 
				for (unsigned int j = 0; j < obj.size(); j++) 
				{
					ReconstructedParticle * reco = dynamic_cast< ReconstructedParticle * >(obj[j]);
					if (weights[j] > maxweight && std::abs(reco->getCharge()) > 0.09) 
					{
						maxweight = weights[j];
						winner = reco;
					}
				}
				if (winner) 
				{
					result->push_back(winner);
				}
			}
		}
		return result;
	}

	vector< Jet * > * JetVertexOperator::TagJets(LCCollection * jetcol, LCCollection * mccol, LCCollection * rel)
	{
		int mcnumber = mccol->getNumberOfElements();
		int jetnumber = jetcol->getNumberOfElements();
		std::cout << "#jets: " << jetnumber << " #vtx: " << mcnumber << '\n';
		vector< Jet * > * result = new vector< Jet * >();
		LCRelationNavigator navigator(rel);
		PIDHandler pidh(jetcol);
		int alid;
		try
		{
			alid = pidh.getAlgorithmID("vtxrec");
		}
		catch(UTIL::UnknownAlgorithm &e)
		{
			std::cout << "No algorithm vtxrec!\n";
			alid = pidh.getAlgorithmID(myAlgorithmName);
		}	
		int taken = -1;
		for (int i = 0; i < mcnumber; i++) 
		{
			Vertex * mcvertex = dynamic_cast< Vertex * >(mccol->getElementAt(i));
			double * positionmc = MathOperator::toDoubleArray(mcvertex->getPosition(),3);
			std::cout << "Vertex tags: " << mcvertex->getParameters()[1] << " & " << mcvertex->getParameters()[2] << '\n';
			if (std::abs(mcvertex->getParameters()[2]) == 2) 
			{
				ReconstructedParticle * mcparticle = mcvertex->getAssociatedParticle();
				//bool found = false;
				float minangle = myAngleCut;
				int winner = -1;
				for (int j = 0; j < jetnumber; j++) 
				{
					ReconstructedParticle * jetpart = dynamic_cast< ReconstructedParticle * >(jetcol->getElementAt(j));
					if (taken == jetpart->id()) 
					{
						std::cout << "Jet skipped by index\n";
						continue;
					}
					float angle = MathOperator::getAngleBtw(mcparticle->getMomentum(), jetpart->getMomentum());
					std::cout << "Angle: " << angle << "\n";
					if (angle < minangle)
					{
						minangle = angle;
						winner = j;
					}//*/
					/*vector< Vertex * > * vertices = convert(navigator.getRelatedToObjects(jetpart));
					const vector< ReconstructedParticle * > components = jetpart->getParticles();
					int nvtx = vertices->size();
					if (nvtx > 0) 
					{
						Vertex * recovertex = vertices->at(0);
						double * positionreco = MathOperator::toDoubleArray(recovertex->getPosition(),3);
						float angle = MathOperator::getAngleBtw(positionreco, positionmc);
						std::cout << "Angle: " << angle << "\n";
						if (angle < minangle)  
						{
							minangle = angle;
							winner = j;
						}
					}//*/
				}
				/*minangle = myAngleCut;
				for (int j = 0; j < jetnumber; j++) 
				{
					if (winner > 0) 
					{
						break;
					}
					ReconstructedParticle * jetpart = dynamic_cast< ReconstructedParticle * >(jetcol->getElementAt(j));
					vector< Vertex * > * vertices = convert(navigator.getRelatedToObjects(jetpart));
					const vector< ReconstructedParticle * > components = jetpart->getParticles();
					int nvtx = vertices->size();
					if (nvtx == 0) 
					{
						for (int k = 0; k < components.size(); k++) 
						{
							ReconstructedParticle * leader = components.at(k);
							float angle = MathOperator::getAngle(mcparticle->getMomentum(), leader->getMomentum());
							if (angle < minangle) 
							{
								minangle = angle;
								winner = j;
							}
							
						}
					}
				}//*/
				if (winner > -1) 
				{
					ReconstructedParticle * jetpart = dynamic_cast< ReconstructedParticle * >(jetcol->getElementAt(winner));
					vector< Vertex * > * vertices = convert(navigator.getRelatedToObjects(jetpart));
					const vector< ReconstructedParticle * > components = jetpart->getParticles();
					int nvtx = vertices->size();
					const ParticleID& pid = pidh.getParticleID(jetpart,alid);
					vector<float> params = pid.getParameters();
					float btag = params[pidh.getParameterIndex(alid,"BTag")];
					float ctag = params[pidh.getParameterIndex(alid,"CTag")];
					taken = jetpart->id();
					std::cout << "Jet energy: " << jetpart->getEnergy() 
						  << " b-tag: " << btag 
						  << " c-tag: " << ctag 
						  << " angle-tag: " << minangle 
						  << " # of vtx: " << nvtx 
						  <<  '\n';
					if (nvtx > 0)
					{
						std::cout << "Algorithm type: " << vertices->at(0)->getAlgorithmType() << "\n";
					}
					Jet * jet = new Jet(btag, ctag, nvtx, jetpart->getMomentum());
					//jet->SetBTag(btag);
					//jet->SetCTag(ctag);
					jet->SetMCPDG(mcvertex->getParameters()[1]);
					jet->SetParticles(components);
					jet->SetRecoVertices(vertices);
					result->push_back(jet);
				}
				else 
				{
					std::cout << "ERROR: Jet not found!!!\n";
				}
			}
		}

		return result;
	}
	bool JetVertexOperator::compareToVertex(Vertex * mcvertex,  vector< Vertex * > * vertices)
	{
		for (unsigned int i = 0; i < vertices->size(); i++) 
		{
			Vertex * recovertex = vertices->at(i);
			float angle = MathOperator::getAngle(mcvertex->getAssociatedParticle()->getMomentum(), recovertex->getAssociatedParticle()->getMomentum());
			if (angle < myAngleCut / 2.0) 
			{
				return true;
			}
		}
		return false;
	}
	void JetVertexOperator::TagVertices(std::vector< Jet * > * jets, EVENT::LCCollection * mccol, vector< Vertex * > & lost)
	{
		int mcnumber = mccol->getNumberOfElements();

		vector<Vertex*> bmcvertices;
		vector<Vertex*> bbarmcvertices;
		for (int i = 0; i < mcnumber; i++) 
		{
			Vertex* vertex = dynamic_cast< Vertex * >(mccol->getElementAt(i));
			if (vertex->getParameters()[1] == 5) 
			{
				bmcvertices.push_back(vertex);
			}
			if (vertex->getParameters()[1] == -5) 
			{
				bbarmcvertices.push_back(vertex);
			}
		}
		bool bexists = false;
		bool bbarexists = false;

		for (unsigned int i = 0; i < jets->size(); i++) 
		{
			Jet * jet = jets->at(i);
			if (jet->GetMCPDG() == 5) 
			{
				bexists = true;
				produceTags(jet, bmcvertices);
			}
			if (jet->GetMCPDG() == -5) 
			{
				bbarexists = true;
				produceTags(jet, bbarmcvertices);
			}
		}

		if (!bexists && bmcvertices.size() > 0) 
		{
			for (unsigned int i = 0; i < bmcvertices.size(); i++) 
			{
				lost.push_back(bmcvertices[i]);
			}
		}
		if (!bbarexists && bbarmcvertices.size() > 0) 
		{
			for (unsigned int i = 0; i < bbarmcvertices.size(); i++)
			{
			        lost.push_back(bbarmcvertices[i]);
			}
		}
	}
	
	vector< ReconstructedParticle * > * JetVertexOperator::GetMissedTracksRel(LCCollection * rel, vector< Jet * > * jets, std::vector< Particle > * converted, LCCollection * out)
	{
		vector< ReconstructedParticle * > * result = new vector< ReconstructedParticle * >();
		std::cout << "Start to extruct missing tracks...\n";
		for (unsigned int i = 0; i < jets->size(); i++) 
		{
			std::cout << "Jet has " << jets->at(i)->GetVertexTags().size() << " tags\n";
			float btag = jets->at(i)->GetBTag();
			if (jets->at(i)->GetVertexTags().size() == 1) 
			{
				VertexTag * tag1 = jets->at(i)->GetVertexTags()[0];
				if (tag1->GetStatus() != EMPTY_TAG) 
				{
					CompareTracks(tag1->__GetMCVertex(), tag1->GetVertex()->getAssociatedParticle()->getParticles(), result,rel, converted, btag);
				}
				continue;
			}
			VertexTag * tag1 = jets->at(i)->GetVertexTags()[0];
			VertexTag * tag2 = jets->at(i)->GetVertexTags()[1];
			if (tag2->GetStatus() == MERGED_TAG && tag1->GetStatus() == MERGED_TAG) 
			{
				std::cout << "We have a merged vertex case with " << tag1->__GetMCVertex()->getAssociatedParticle()->getParticles().size() + tag2->__GetMCVertex()->getAssociatedParticle()->getParticles().size() <<" mcparticles\n";
				CompareTracks(tag1->__GetMCVertex(), tag1->GetVertex()->getAssociatedParticle()->getParticles(), result,rel, converted,btag);
				CompareTracks(tag2->__GetMCVertex(), tag2->GetVertex()->getAssociatedParticle()->getParticles(), result,rel, converted,btag);
			}
			if (tag2->GetStatus() == PRECISE_TAG && tag1->GetStatus() == PRECISE_TAG) 
			{
				const vector< ReconstructedParticle * > fromtag1reco = tag1->GetVertex()->getAssociatedParticle()->getParticles();
				const vector< ReconstructedParticle * > fromtag2reco = tag2->GetVertex()->getAssociatedParticle()->getParticles();
				vector< ReconstructedParticle * > recounited;
				recounited.reserve(fromtag1reco.size() + fromtag2reco.size());
				recounited.insert(recounited.end(),fromtag1reco.begin(),fromtag1reco.end());
				recounited.insert(recounited.end(),fromtag2reco.begin(),fromtag2reco.end());
				std::cout << "We have a 2 recovertices case with " <<  recounited.size() << " recoparticles\n";
				CompareTracks(tag1->__GetMCVertex(), recounited, result,rel, converted, btag);
				CompareTracks(tag2->__GetMCVertex(), recounited, result,rel, converted, btag);
			}
			if (tag1->GetStatus() == PRECISE_TAG && tag2->GetStatus() == EMPTY_TAG) 
			{
				
				std::cout << "We have one vertex case with " << tag1->__GetMCVertex()->getAssociatedParticle()->getParticles().size() << " mc particles & " << tag1->__GetMCVertex()->getAssociatedParticle()->getParticles().size() << " recoparticles\n";
				CompareTracks(tag1->__GetMCVertex(), tag1->GetVertex()->getAssociatedParticle()->getParticles(), result, rel, converted,btag);
			}
			if (tag1->GetStatus() == EMPTY_TAG && tag2->GetStatus() == PRECISE_TAG) 
			{
				std::cout << "We have one vertex case with " << tag2->__GetMCVertex()->getAssociatedParticle()->getParticles().size() << " mc particles & " << tag2->__GetMCVertex()->getAssociatedParticle()->getParticles().size() << " recoparticles\n";
				CompareTracks(tag2->__GetMCVertex(), tag2->GetVertex()->getAssociatedParticle()->getParticles(), result, rel, converted,btag);
			}
		}
		return result;
	}

	void JetVertexOperator::CompareTracks(Vertex * mcvertex, const vector< ReconstructedParticle * > & recotracks, vector< ReconstructedParticle * > * missedTotal, LCCollection * rel, vector< Particle > * convertedTotal, float btag)
	{
		const vector< ReconstructedParticle * > mctracks = mcvertex->getAssociatedParticle()->getParticles();
		vector< ReconstructedParticle * > recoPFOtracks = mapToPFO(recotracks);
		vector< MCParticle * > recoRelatedPFOtracks = ParticleOperator::GetMCParticlesRel(recoPFOtracks, rel, myTrackRel);

		//vector< Particle > * result = new vector<Particle >();
		vector< ReconstructedParticle * > * missed = CompareTracksRel(recoRelatedPFOtracks, mctracks);
		vector< MCParticle * > egmissed = mapToProngs(*missed);
		
		missedTotal->reserve(missedTotal->size() +  missed->size());
		missedTotal->insert(missedTotal->end(),missed->begin(),missed->end());
		std::cout << "Number of missed tracks: " << missedTotal->size()  << '\n';
		if (!convertedTotal) 
		{
			return;
		}
		for (int i = 0; i < missed->size(); i++) 
		{
			ReconstructedParticle * particle = missed->at(i);

			MCParticle * missprong = egmissed.at(i);
			ReconstructedParticle * matchedmissed = ParticleOperator::GetRecoParticle(missprong, rel);

			vector< float > direction = MathOperator::getDirection(particle->getMomentum());
			double * sec = MathOperator::toDoubleArray(mcvertex->getPosition(), 3);
			vector< float > angles = MathOperator::getAngles(direction);
			float offset = MathOperator::getDistanceTo(ip, direction, sec);
			Particle missedParticle;
			missedParticle.SetJetBtag(btag);
			missedParticle.SetOffset(offset);
			missedParticle.SetTheta(angles[1]);
			missedParticle.SetMomentum(particle->getMomentum());
			missedParticle.SetMass(particle->getMass());
			missedParticle.SetVertex(mcvertex);
			missedParticle.SetGeneration(mcvertex->getParameters()[2]);
			missedParticle.SetInteracted(missprong->isDecayedInCalorimeter());
			missedParticle.SetMCParticle(missprong);
			int * hits = new int[4];
			if (matchedmissed) 
			{
				missedParticle.SetVertexAngle(getVertexAngle(matchedmissed, mcvertex));
				missedParticle.SetChi2(matchedmissed->getTracks()[0]->getChi2() / (float) matchedmissed->getTracks()[0]->getNdf());
				missedParticle.SetIsReco(1);
				hits[0] = matchedmissed->getTracks()[0]->getSubdetectorHitNumbers()[0];
				hits[1] = matchedmissed->getTracks()[0]->getSubdetectorHitNumbers()[4];
				hits[2] = matchedmissed->getTracks()[0]->getSubdetectorHitNumbers()[6];
				hits[3] = matchedmissed->getTracks()[0]->getSubdetectorHitNumbers()[2];
				missedParticle.SetTruthAngle(MathOperator::getAngleBtw(missprong->getMomentum(), matchedmissed->getMomentum()));
				missedParticle.SetHits(hits);
				missedParticle.SetRecoParticle(matchedmissed);
			}
			else 
			{
				missedParticle.SetChi2(-1.0);
				missedParticle.SetIsReco(0);
			}
			convertedTotal->push_back(missedParticle);
			std::cout << "INFO: MISSED TRACK OFFSET:  " << offset
				  << " VXD: " << hits[0]
				  << " MOMENTUM: " << MathOperator::getModule(missedParticle.GetMomentum())
				  << " COSTHETA: " << std::cos(missedParticle.GetTheta())
				  << '\n';
		}
	}

	std::vector< EVENT::MCParticle * > * JetVertexOperator::GetMissedTracks(EVENT::LCCollection * prongs, EVENT::LCCollection * rel,  EVENT::LCCollection * out)
	{
		vector< MCParticle * > * result = new vector< MCParticle * >();
		std::cout << "Start to extruct all missing tracks...\n";
		int prongnumber = prongs->getNumberOfElements();
		LCRelationNavigator navigator(rel);
		vector< int > trackIDs;
		prongs->getParameters().getIntVals("trackIDs",trackIDs);
		for (int i = 0; i < prongnumber; i++) 
		{
			MCParticle * particle =  dynamic_cast< MCParticle * >(prongs->getElementAt(i));
			//std::cout << "Processing particle with id " << trackIDs[i] << " and momentum " << MathOperator::getModule(particle->getMomentum()) << " GeV:\n";
			vector< LCObject * > obj = navigator.getRelatedFromObjects(particle);
			if (obj.size() < 1) 
			{
				std::cout << "INFO: Lost truthlink of particle with momentum of " << MathOperator::getModule(particle->getMomentum()) << " GeV\n";
				if (out) 
				{
					out->addElement(particle);
				}
				result->push_back(particle);
				continue;
			}
			if (obj.size() > 0) 
			{
				bool matched = false;
				int winner = -1;
				float maxweight = 0.0;
				vector< float > weights = navigator.getRelatedFromWeights (particle); 
				for (unsigned int j = 0; j < obj.size(); j++) 
				{
					ReconstructedParticle * reco = dynamic_cast< ReconstructedParticle * >(obj[j]);
					if (weights[j] > maxweight && std::abs(reco->getCharge()) > 0.09) 
					{
						maxweight = weights[j];
						winner = j;
					}
				}
				ReconstructedParticle * reco = dynamic_cast< ReconstructedParticle * >(obj[(winner == -1)? 0:winner]);
				if (std::abs(reco->getCharge()) < 0.09) 
				{
					std::cout << "FATAL: " << weights[(winner == -1)? 0:winner] 
						<< " truthlink for " << MathOperator::getModule(particle->getMomentum())
						<< " particle points to a neutral "
						<< MathOperator::getModule(reco->getMomentum()) 
						<< " recoparticle!  Charges: "
						<< particle->getCharge() << " and " 
						<< reco->getCharge() << "\n";
					if (obj.size() == 1) 
					{
						std::cout << "FATAL: and... It is the only truthlink!\n";
					}
					if (out)
					{
					        out->addElement(particle);
					}
					result->push_back(particle);
					continue;
				}
				else 
				{
					/*std::cout << "INFO: " << weights[winner] 
						<< " truthlink for " << MathOperator::getModule(particle->getMomentum())
						<< " particle points to a "
						<< MathOperator::getModule(reco->getMomentum()) 
						<< " recoparticle.  Charges: "
						<< particle->getCharge() << " and " 
						<< reco->getCharge() << "\n";*/
					
				}
				if (!ParticleOperator::CheckForVertex(reco) || winner == -1) 
				{
					std::cout << "INFO: Lost particle with momentum of " 
					<< MathOperator::getModule(particle->getMomentum()) 
					<< " GeV! It points to " 
					<< MathOperator::getModule(reco->getMomentum()) 
					<< " GeV recoparticle.\n";
					if (out)
					{
					        out->addElement(particle);
					}
					result->push_back(particle);

				}
			}
		}
		return result;
	}
	vector< ReconstructedParticle * > * JetVertexOperator::GetMissedTracks(std::vector< Jet * > * jets, std::vector< Particle > * converted)
	{
		vector< ReconstructedParticle * > * result = new vector< ReconstructedParticle * >();
		std::cout << "Start to extruct missing tracks...\n";
		for (unsigned int i = 0; i < jets->size(); i++) 
		{
			std::cout << "Jet has " << jets->at(i)->GetVertexTags().size() << " tags\n";
			VertexTag * tag1 = jets->at(i)->GetVertexTags()[0];
			VertexTag * tag2 = jets->at(i)->GetVertexTags()[1];
			if (tag2->GetStatus() == MERGED_TAG && tag1->GetStatus() == MERGED_TAG) 
			{
				std::cout << "We have a merged vertex case with " << tag1->__GetMCVertex()->getAssociatedParticle()->getParticles().size() + tag2->__GetMCVertex()->getAssociatedParticle()->getParticles().size() <<" mcparticles\n";
				CompareTracks(tag1->__GetMCVertex(), tag1->GetVertex()->getAssociatedParticle()->getParticles(), result, converted);
				CompareTracks(tag2->__GetMCVertex(), tag2->GetVertex()->getAssociatedParticle()->getParticles(), result, converted);
			}
			if (tag2->GetStatus() == PRECISE_TAG && tag1->GetStatus() == PRECISE_TAG) 
			{
				const vector< ReconstructedParticle * > fromtag1reco = tag1->GetVertex()->getAssociatedParticle()->getParticles();
				const vector< ReconstructedParticle * > fromtag2reco = tag2->GetVertex()->getAssociatedParticle()->getParticles();
				vector< ReconstructedParticle * > recounited;
				recounited.reserve(fromtag1reco.size() + fromtag2reco.size());
				recounited.insert(recounited.end(),fromtag1reco.begin(),fromtag1reco.end());
				recounited.insert(recounited.end(),fromtag2reco.begin(),fromtag2reco.end());
				std::cout << "We have a 2 recovertices case with " <<  recounited.size() << " recoparticles\n";
				CompareTracks(tag1->__GetMCVertex(), recounited, result, converted);
				CompareTracks(tag2->__GetMCVertex(), recounited, result, converted);
			}
			if (tag1->GetStatus() == PRECISE_TAG && tag2->GetStatus() == EMPTY_TAG) 
			{
				
				std::cout << "We have one vertex case with " << tag1->__GetMCVertex()->getAssociatedParticle()->getParticles().size() << " mc particles & " << tag1->__GetMCVertex()->getAssociatedParticle()->getParticles().size() << " recoparticles\n";
				CompareTracks(tag1->__GetMCVertex(), tag1->GetVertex()->getAssociatedParticle()->getParticles(), result, converted);
			}
			if (tag1->GetStatus() == EMPTY_TAG && tag2->GetStatus() == PRECISE_TAG) 
			{
				std::cout << "We have one vertex case with " << tag2->__GetMCVertex()->getAssociatedParticle()->getParticles().size() << " mc particles & " << tag2->__GetMCVertex()->getAssociatedParticle()->getParticles().size() << " recoparticles\n";
				CompareTracks(tag2->__GetMCVertex(), tag2->GetVertex()->getAssociatedParticle()->getParticles(), result, converted);
			}
		}
		return result;
	}
	void JetVertexOperator::CompareTracks(EVENT::Vertex * mcvertex, const std::vector< EVENT::ReconstructedParticle * > & recotracks, std::vector< EVENT::ReconstructedParticle * > * missedTotal, std::vector< Particle > * convertedTotal)
	{
		const vector< ReconstructedParticle * > mctracks = mcvertex->getAssociatedParticle()->getParticles();
		//vector< Particle > * result = new vector<Particle >();
		vector< ReconstructedParticle * > * missed = CompareTracks(recotracks, mctracks);
		missedTotal->reserve(missedTotal->size() +  missed->size());
		missedTotal->insert(missedTotal->end(),missed->begin(),missed->end());
		std::cout << "Number of missed tracks: " << missedTotal->size()  << '\n';
		if (!convertedTotal) 
		{
			return;
		}
		for (int i = 0; i < missed->size(); i++) 
		{
			ReconstructedParticle * particle = missed->at(i);
			vector< float > direction = MathOperator::getDirection(particle->getMomentum());
			double * sec = MathOperator::toDoubleArray(mcvertex->getPosition(), 3);
			vector< float > angles = MathOperator::getAngles(direction);
			float offset = MathOperator::getDistanceTo(ip, direction, sec);
			Particle missedParticle;
			missedParticle.SetOffset(offset);
			missedParticle.SetTheta(angles[1]);
			missedParticle.SetMomentum(particle->getMomentum());
			missedParticle.SetMass(particle->getMass());
			missedParticle.SetVertex(mcvertex);
			convertedTotal->push_back(missedParticle);
			std::cout << "INFO: MISSED TRACK OFFSET:  " << offset
				  << " MOMENTUM: " << MathOperator::getModule(missedParticle.GetMomentum())
				  << " THETA: " << missedParticle.GetTheta()
				  << '\n';
		}
	}

	vector< ReconstructedParticle * > * JetVertexOperator::CompareTracks(const vector< ReconstructedParticle * > & recotracks, const vector< ReconstructedParticle * > & mctracks)
	{
		vector< ReconstructedParticle * > * result = new vector< ReconstructedParticle * >();
		for (unsigned int i = 0; i < mctracks.size(); i++) 
		{
			ReconstructedParticle * mcparticle = mctracks[i];
			bool found = false;
			for (unsigned int j = 0; j < recotracks.size(); j++) 
			{
				ReconstructedParticle * recoparticle = recotracks[j];
				if (ParticleOperator::CompareParticles(mcparticle,recoparticle)) 
				{
					found = true;
					break;
				}
			}
			if (!found) 
			{
				result->push_back(mcparticle);
			}
		}
		return result;
	}
	vector< ReconstructedParticle * > * JetVertexOperator::CompareTracksRel(const vector< MCParticle * > & recotracks, const vector< ReconstructedParticle * > & mctracks)
	{
		vector< ReconstructedParticle * > * result = new vector< ReconstructedParticle * >();
		for (unsigned int i = 0; i < mctracks.size(); i++) 
		{
			ReconstructedParticle * mcparticle = mctracks[i];
			bool found = false;
			for (unsigned int j = 0; j < recotracks.size(); j++) 
			{
				MCParticle * recoparticle = recotracks[j];
				if (ParticleOperator::CompareParticles(recoparticle,mcparticle)) 
				{
					found = true;
					break;
				}
			}
			if (!found) 
			{
				result->push_back(mcparticle);
			}
		}
		return result;
	}
	
	vector< VertexTag * > * JetVertexOperator::tagOneVertex(vector< Vertex * > & mcvertices, Vertex * recovertex)
	{
		std::cout << "\tFound only one recovertex!\n";
		float distance1 = MathOperator::getDistance(mcvertices[0]->getPosition(), recovertex->getPosition());
		float distance2 = MathOperator::getDistance(mcvertices[1]->getPosition(), recovertex->getPosition());
		int reconumber = recovertex->getAssociatedParticle()->getParticles().size();
		int mcnumber1 = mcvertices[0]->getAssociatedParticle()->getParticles().size();
		int mcnumber2 = mcvertices[1]->getAssociatedParticle()->getParticles().size();
		float distance = (distance1 + distance2);
		float ratio = distance1 / distance;
		std::cout << "\tThe mc distances are: " << distance1 << " & " << distance2
			  << " the mc numbers are: " << mcnumber1<< " & " << mcnumber2 
			  << " the ratio of both is: " << ratio << '\n'
			  << "\tThe reco number is: " << reconumber << '\n';
		VertexTag * tag = NULL;
		VertexTag * other = NULL;
		if ((distance1 < distance2 && ratio < 0.1 && reconumber <= mcnumber1) && distance > 1.0)
		{
			tag = new VertexTag(mcvertices[0]);
			tag->SetRecoVertex(recovertex);
			tag->SetStatus(PRECISE_TAG);
			other =  new VertexTag(mcvertices[1]);
			other->SetStatus(EMPTY_TAG);
			std::cout << "\tThe vertex with " << mcnumber1 << " tracks is tagged, other one is empty\n";
		}
		if ((distance2 < distance1 && ratio > 0.9 && reconumber <= mcnumber2) && distance > 1.0)
		{
			tag = new VertexTag(mcvertices[1]);
			tag->SetRecoVertex(recovertex);
			tag->SetStatus(PRECISE_TAG);
			other =  new VertexTag(mcvertices[0]);
			other->SetStatus(EMPTY_TAG);
			std::cout << "\tThe vertex with " << mcnumber2 << " tracks is tagged, other one is empty\n";
		}
		if (!tag) 
		{
			std::cout << "\tThe difficult case: vertex marked as merged!\n";
			tag = new VertexTag(mcvertices[0]);
			tag->SetRecoVertex(recovertex);
			other =  new VertexTag(mcvertices[1]);
			other->SetRecoVertex(recovertex);
			other->SetStatus(MERGED_TAG);
			tag->SetStatus(MERGED_TAG);
		}
		if (!tag) 
		{
			std::cout << "\tERROR: Vertices remain untagged!\n";
			return NULL;
		}
		vector< VertexTag * > * result = new vector< VertexTag * >();
		result->push_back(tag);
		result->push_back(other);
		return result;
	}
	
	void JetVertexOperator::produceTags(Jet * jet, std::vector< Vertex * > & mcvertices, std::vector< VertexTag * > * tags)
	{
		vector< Vertex * > * recovertices = jet->GetRecoVertices();
		if (mcvertices.size() < 1) 
		{
			std::cout << "MC vertices are uninitialized!\n";
			return;
		}
		int taken = -1;
		std::cout << "Begin tagging:\n";
		if (mcvertices.size() == 1) 
		{
			Vertex * mcvertex = mcvertices[0];
			std::cout << "\tOne mc vertex mode:\n";
			if (recovertices->size() == 0) 
			{
				VertexTag * tag = new VertexTag(mcvertex);
				tag->SetStatus(EMPTY_TAG);
				jet->AddVertexTag(tag);
			}
			for (int i = 0; i < recovertices->size(); i++) 
			{
				VertexTag * tag = new VertexTag(mcvertex);
				tag->SetRecoVertex(recovertices->at(i));
				tag->SetStatus(MERGED_TAG);
				jet->AddVertexTag(tag);
			}
			return;
		}
		if (recovertices->size() == 1) 
		{
			Vertex * recovertex = recovertices->at(0);
			vector< VertexTag * > * ttags = tagOneVertex(mcvertices,recovertex);
			//if (tags) 
			{
				jet->AddVertexTag(ttags->at(0));
				jet->AddVertexTag(ttags->at(1));
			}
			return;
		}
		for (unsigned int i = 0; i < mcvertices.size(); i++) 
		{
			Vertex * mcvertex = mcvertices[i];
			VertexTag * tag = new VertexTag(mcvertex);
			if (recovertices->size() == 0) 
			{
				std::cout << "\tTag is empty - no recovertices\n";
				tag->SetStatus(EMPTY_TAG);
				jet->AddVertexTag(tag);
			}
			if (recovertices->size() == 2 && taken < 0) 
			{
				float distance1 = MathOperator::getDistance(mcvertex->getPosition(), recovertices->at(0)->getPosition());
				float distance2 = MathOperator::getDistance(mcvertices[1]->getPosition(), recovertices->at(1)->getPosition());
				if (distance1 < distance2) 
				{
					tag->SetRecoVertex(recovertices->at(1)); // CHANGED
					taken = 1; // CHANGED
				}
				else 
				{
					tag->SetRecoVertex(recovertices->at(0));
					taken = 0;
				}
				std::cout << "\tOne of two vertices tagged with " << distance1 << " vs " << distance2 <<" \n";
				tag->SetStatus(PRECISE_TAG);//PRECISE_TAG
				jet->AddVertexTag(tag);
				continue;
			}
			if (recovertices->size() == 2 && taken > -1) 
			{
				tag->SetRecoVertex((taken)?recovertices->at(0):recovertices->at(1));
				std::cout << "\tAnother one tagged\n";
				tag->SetStatus(PRECISE_TAG);
				jet->AddVertexTag(tag);
			}
			if (tags) 
			{
				tags->push_back(tag);
			}
		}
	}
	std::vector< Vertex * > * JetVertexOperator::convert(const std::vector< LCObject * > & objs)
	{
		std::vector< Vertex * > * result = new std::vector< Vertex * >();
		for (int i = 0; i < objs.size(); i++) 
		{
			result->push_back(dynamic_cast< Vertex * >(objs[i]));
		}
		return result;
	}
	
	vector< ReconstructedParticle * > JetVertexOperator::mapToPFO(const vector< ReconstructedParticle * > & secondaries)
	{
		int pfonumber = myPFO->getNumberOfElements();
		vector< ReconstructedParticle * > result;
		/*for (int i = 0; i < pfonumber; i++) 
		{
			ReconstructedParticle * recopfo = dynamic_cast< ReconstructedParticle * >(myPFO->getElementAt(i));
			if (std::abs(recopfo->getCharge()) < 0.1 ) 
			{
				continue;
			}
			//std::cout << "Particle mapping begin:\n";
			bool mapped = false;
			for (int j = 0; j < secondaries.size(); j++) 
			{
				if (ParticleOperator::CompareParticles(recopfo, secondaries[j])) 
				{
					result.push_back(recopfo);
					mapped = true;
					//std::cout << "Particle mapped!\n";
					break;
				}
			}
			if (!mapped) 
			{
				std::cout << "ERROR: Particle is NOT mapped!\n";
			}
		}*/
		for (int j = 0; j < secondaries.size(); j++) 
		{
			bool mapped = false;
			for (int i = 0; i < pfonumber; i++) 
			{
				ReconstructedParticle * recopfo = dynamic_cast< ReconstructedParticle * >(myPFO->getElementAt(i));
				if (std::abs(recopfo->getCharge()) < 0.1 )
				{
				        continue;
				}
				if (ParticleOperator::CompareParticles(recopfo, secondaries[j]))
				{
					result.push_back(recopfo);
					mapped = true;
					//std::cout << "Particle mapped!\n";
					break;
				}

			}
			if (!mapped) 
			{
			        std::cout << "ERROR: Particle is NOT mapped!\n";
				result.push_back(secondaries[j]);
			}
		}
		return result;
	}
	vector< MCParticle * > JetVertexOperator::mapToProngs(const vector< ReconstructedParticle * > & secondaries)
	{
		int number = myEGProngs->getNumberOfElements();
		vector< MCParticle * > result;
		for (int i = 0; i < number; i++) 
		{
			MCParticle * prong = dynamic_cast< MCParticle * >(myEGProngs->getElementAt(i));
			for (int j = 0; j < secondaries.size(); j++) 
			{
				if (ParticleOperator::CompareParticles(prong, secondaries[j])) 
				{
					result.push_back(prong);
					break;
				}
			}
		}
		return result;
	}
	float JetVertexOperator::getVertexAngle(ReconstructedParticle * particle, Vertex * secvertex)
	{
		/*double * start = new double[3];
		Track * track = NULL; //particle->getTracks()[0];
		for (int i = 0; i < particle->getTracks()[0]->getTracks().size(); i++) 
		{
			if (particle->getTracks()[0]->getTracks()[i]->getSubdetectorHitNumbers()[0] > 0) 
			{
				std::cout << " i: " << i << " hits: " << particle->getTracks()[0]->getTracks()[i]->getSubdetectorHitNumbers()[0] << " size: " << particle->getTracks()[0]->getTracks().size() <<  "\n";
				track = particle->getTracks()[0]->getTracks()[i];
				particle = ParticleOperator::ReconstructParticle(track);
			}
		}
		if (!track) 
		{
			track = particle->getTracks()[0];
		}
		float d0 = track->getD0();
		float z0 = track->getZ0();
		float phi0 = track->getPhi();
		start[0] =  - d0 * std::sin(phi0);
		start[1] =  d0 * std::cos(phi0);
		start[2] = z0;
		*/
		if (!particle || particle->getTracks().size() < 1) 
		{
			return -1.;
		}
		TrackOperator opera; 
		double * start = opera.GetStartPoint(particle);
		double * secondaryPosition = MathOperator::toDoubleArray(secvertex->getPosition(),3);
		vector<float> diff = MathOperator::getDirection(secondaryPosition, start);
		double diif[3];
		for (int m = 0; m < 3; m++) 
		{
			diif[m] = diff[m];
		}
		float angle = MathOperator::getAngle(diif, particle->getMomentum());
		return angle;
	}



}

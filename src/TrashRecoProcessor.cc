#include "TrashRecoProcessor.hh"
using std::vector;
using EVENT::Vertex;
using IMPL::VertexImpl;
using TTbarAnalysis::MathOperator;

namespace TTbarAnalysis 
{
	TrashRecoProcessor aTrashRecoProcessor ;
	TrashRecoProcessor::TrashRecoProcessor() : Processor("TrashRecoProcessor") {
	
	    // modify processor description
	    _description = "TrashRecoProcessor does whatever it does ..." ;
	
	
	    // register steering parameters: name, description, class-variable, default value
	    registerInputCollection( LCIO::VERTEX,
	            "SecondaryCollectionName" , 
	            "Name of the BuildUpVertex collection"  ,
	            _colSecName ,
	            std::string("BuildUpVertex")
	    );
	    registerInputCollection( LCIO::VERTEX,
	            "PrimaryCollectionName" , 
	            "Name of the PrimaryVertex collection"  ,
	            _colPriName ,
	            std::string("PrimaryVertex")
	
	    );
	    registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
        	    "OutputMissedName" , 
	            "Name of the missed collection"  ,
        	    _colMissName ,
           	 std::string("MissedParticles")
	    );
	    registerOutputCollection( LCIO::VERTEX,
        	    "OutputTaggedName" , 
	            "Name of the output vertex collection"  ,
        	    _colTagName ,
           	 std::string("TaggedVertices")
	    );
	    registerInputCollection( LCIO::VERTEX,
	    	"MCVertexCollectionName",
		"Name of the MCVertex collection",
		_colMCName,
	    	std::string("MCVertex")
	    );
	    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
	    	"JetCollectionName",
		"Name of the Jet collection",
		_colJetsName,
	    	std::string("FinalJets")
	    );
	    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
	    	"JetRelCollectionName",
		"Name of the Jet relation collection",
		_colJetsRelName,
	    	std::string("FinalJets_rel")
	    );
	    _angleAcceptance = 0.2;
	    registerProcessorParameter("angleAcceptance",
	    	"Angle cut for truth tagging verticies by direction",
		_angleAcceptance,
		_angleAcceptance);
	    _handleJets = 0;
	    registerProcessorParameter("handleJets",
	    	"Handle jets for b-tagging",
		_handleJets,
		_handleJets);
	    registerInputCollection( LCIO::MCPARTICLE,
	    	"QuarksCollectionName",
		"Name of the quarks collection",
		_colquarkName,
	    	std::string("MCbquarks")
	    );

	}
	
	
	
	void TrashRecoProcessor::init() 
	{ 
		streamlog_out(DEBUG) << "   init called  " << std::endl ;
		printParameters() ;
		_nRun = 0 ;
		_nEvt = 0 ;
		_hfilename = "TrashRecoTest.root";
		_hfile = new TFile( _hfilename.c_str(), "RECREATE", _hfilename.c_str() ) ;
		_hTree = new TTree( "Stats", "My vertex tree!" );
		_hTree->Branch("numberOfTagged", &_numberOfTagged, "numberOfTagged/I");
		_hTree->Branch("numberOfTotal", &_numberOfTotal, "numberOfTotal/I");
		_hTree->Branch("numberOfTernary", &_numberOfTernary, "numberOfTernary/I");
		_hTree->Branch("numberOfSecondary", &_numberOfSecondary, "numberOfSecondary/I");
		_hTree->Branch("numberOfUnknown", &_numberOfUnknown, "numberOfUnknown/I");
		_hTree->Branch("bexists", &_bexists, "bexists/I");
		_hTree->Branch("bbarexists", &_bbarexists, "bbarexists/I");
		_hTree->Branch("bnumber", &_bnumber1, "bnumber1/I");
		_hTree->Branch("bbarnumber", &_bbarnumber1, "bbarnumber1/I");
		_hTree->Branch("bcharge", &_bcharge, "bcharge/I");
		_hTree->Branch("bbarcharge", &_bbarcharge, "bbarcharge/I");
		_hTree->Branch("bteta", &_bteta, "bteta/F");
		_hTree->Branch("bbarteta", &_bbarteta, "bbarteta/F");
		_hTree->Branch("bptmiss", &_bptmiss, "bptmiss/F");
		_hTree->Branch("bbarptmiss", &_bbarptmiss, "bbarptmiss/F");
		_hTree->Branch("bmomentum", &_bmomentum, "bmomentum/F");
		_hTree->Branch("bbarmomentum", &_bbarmomentum, "bbarmomentum/F");
		_hTree->Branch("bbarIPdistance", &_bbarIPdistance, "bbarIPdistance/F");
		_hTree->Branch("bIPdistance", &_bIPdistance, "bIPdistance/F");
		_hTree->Branch("bbarprobmean", &_bbarprobmean, "bbarprobmean/F");
		_hTree->Branch("bprobmean", &_bprobmean, "bprobmean/F");
		_hTree->Branch("bbarchimean", &_bbarchimean, "bbarchimean/F");
		_hTree->Branch("bchimean", &_bchimean, "bchimean/F");
		//_hTree->Branch("numberOfDistances", &_numberOfDistances, "numberOfDistances/I");
		//_hTree->Branch("distances", _distances, "distances[numberOfDistances]/F");
		_hTaggedTree = new TTree( "TaggedVertices", "My vertex tree!" );
		_hTaggedTree->Branch("numberOfTagged", &_numberOfTagged, "numberOfTagged/I");
		_hTaggedTree->Branch("distance", _distanceFromIP, "distance[numberOfTagged]/F");
		_hTaggedTree->Branch("mindistance", _mindistance, "mindistance[numberOfTagged]/F");
		_hTaggedTree->Branch("angle", _angle, "angle[numberOfTagged]/F");
		_hTaggedTree->Branch("probability", _probability, "probability[numberOfTagged]/F");
		_hTaggedTree->Branch("chi2", _chi2, "chi2[numberOfTagged]/F");
		_hTaggedTree->Branch("coordinates", _coordinates, "coordinates[numberOfTagged][3]/F");
		_hTaggedTree->Branch("PDG", _PDG, "PDG[numberOfTagged]/I");
		_hTaggedTree->Branch("generation", _generation, "generation[numberOfTagged]/I");
		_hTaggedTree->Branch("charge", _charge, "charge[numberOfTagged]/I");
		_hTaggedTree->Branch("numberOfParticles", _numberOfParticles, "numberOfParticles[numberOfTagged]/I");
		_hTaggedTree->Branch("energyOfParticles", _energyOfParticles, "energyOfParticles[numberOfTagged][15]/F");
		_hTaggedTree->Branch("momentumOfParticles", _momentumOfParticles, "momentumOfParticles[numberOfTagged][15]/F");
		_hTaggedTree->Branch("massOfParticles", _massOfParticles, "massOfParticles[numberOfTagged][15]/F");
		
		_hUntaggedTree = new TTree( "UntaggedVertices", "My vertex tree!" );
		_hUntaggedTree->Branch("numberOfUnknown", &_numberOfUnknown, "numberOfUnknown/I");
		_hUntaggedTree->Branch("distance", _distanceFromIP, "distance[numberOfUnknown]/F");
		_hUntaggedTree->Branch("probability", _probability, "probability[numberOfUnknown]/F");
		_hUntaggedTree->Branch("chi2", _chi2, "chi2[numberOfUnknown]/F");
		_hUntaggedTree->Branch("charge", _charge, "charge[numberOfUnknown]/I");
		_hUntaggedTree->Branch("coordinates", _coordinates, "coordinates[numberOfUnknown][3]/F");
		_hUntaggedTree->Branch("numberOfParticles", _numberOfParticles, "numberOfParticles[numberOfUnknown]/I");
		_hUntaggedTree->Branch("energyOfParticles", _energyOfParticles, "energyOfParticles[numberOfUnknown][15]/F");
		_hUntaggedTree->Branch("momentumOfParticles", _momentumOfParticles, "momentumOfParticles[numberOfUnknown][15]/F");
		_hUntaggedTree->Branch("massOfParticles", _massOfParticles, "massOfParticles[numberOfUnknown][15]/F");
		//_hTaggedTree->Branch("particles", _particles);
		_hJetTree = new TTree( "Jets", "My vertex tree!" );
		_hJetTree->Branch("numberOfJets", &_numberOfJets, "numberOfJets/I");
		_hJetTree->Branch("numberOfVertices", _nvertices, "numberOfVertices[numberOfJets]/I");
		_hJetTree->Branch("btags", _btags, "btags[numberOfJets]/F");
		_hJetTree->Branch("ctags", _ctags, "ctags[numberOfJets]/F");
		_hJetTree->Branch("mcpdg", _mcpdg, "mcpdg[numberOfJets]/I");
		
		_hMissedTree = new TTree( "MissedTracks", "My vertex tree!" );
		_hMissedTree->Branch("numberOfMissed", &_numberOfMissed, "numberOfMissed/I");
		_hMissedTree->Branch("offsetMissed", _offsetMissed, "offsetMissed[numberOfMissed]/F");
		_hMissedTree->Branch("momentumMissed", _momentumMissed, "momentumMissed[numberOfMissed]/F");
		_hMissedTree->Branch("thetaMissed", _thetaMissed, "thetaMissed[numberOfMissed]/F");
		
	
	}
	
	
	void TrashRecoProcessor::processRunHeader( LCRunHeader* run) 
	{ 
		_nRun++ ;
	} 
	
	Vertex * TrashRecoProcessor::FindPrimaryVertex(const LCCollection * col)
	{
		int number = col->getNumberOfElements();
		if (number > 1) 
		{
			std::cout << "WARNING: Pile-up found!\n";
		}
		return dynamic_cast< Vertex * >( col->getElementAt(0) ) ;
	}
	void TrashRecoProcessor::PrintParticle(ReconstructedParticle * particle)
	{
		if (!particle) 
		{
			return;
		}
		std::cout << std::fixed << std::setw( 6 ) << std::setprecision( 3 ) << std::setfill( ' ' );
		int id = 0;
		if (particle->getParticleIDUsed()) 
		{
			std::cout << "Type " << particle->getParticleIDUsed()->getType() << '\n';
			id = particle->getParticleIDs()[0]->getPDG(); 
		}
		std::cout<<"|"<< id <<"\t\t|"<<particle->getMass()<<"\t\t|"<<particle->getCharge()  <<"\t\t|"<<particle->getEnergy() <<"\t\t|\n";
	
	}
	
	void TrashRecoProcessor::processEvent( LCEvent * evt ) 
	{ 
		try
		{
			LCCollection* col = evt->getCollection( _colSecName );
			LCCollection* mc = evt->getCollection( _colMCName );
			int number = col->getNumberOfElements();
			_numberOfTotal = number;
			//std::cout << "Event: " << _nEvt << '\n';
			_nEvt ++ ;
			_primary = FindPrimaryVertex( evt->getCollection( _colPriName ));
			VertexRecoOperator reco(_angleAcceptance, _primary);
			if (_handleJets) 
			{
				try{
					LCCollection* jets = evt->getCollection( _colJetsName );
					LCCollection* rel = evt->getCollection( _colJetsRelName );
					JetOperator jetOperator(0.1 , "lcfiplus");
					//jetOperator.GetBtags(jets);
					vector< Jet * > * btagjets = jetOperator.GetJets(jets, rel);
					LCCollection* quarks = evt->getCollection( _colquarkName );
					jetOperator.CompareDirection(btagjets, quarks);
					Write(btagjets);
					_hJetTree->Fill();
				}
				catch( DataNotAvailableException &e)
				{
					std::cout << "Jets collections are not available!\n";
				}
			}
			//vector< VertexTag * > * tagged = reco.Compare(col, mc);
			vector< VertexTag * > * tagged = reco.CompareDirection(col, mc);
			vector< ReconstructedParticle * > * recomissed = new vector< ReconstructedParticle * >();
			vector< Particle > * missed = reco.CompareTracks(col, mc, recomissed);
			int reconot = 0;
			int mcnot = 0;
			int missednot = missed->size();
			for (int i = 0; i < tagged->size(); i++) 
			{
				reconot += tagged->at(i)->GetTrackNumber();
				mcnot += tagged->at(i)->__GetMCTrackNumber();
			}
			_numberOfTagged = tagged->size();
			_numberOfUnknown = number - _numberOfTagged;
			std::cout << "We have " << number << " vertices total\n";
			WriteTagged(tagged);
			WriteTaggedCollection(evt, tagged);
			Write(evt, missed, recomissed);
			_hMissedTree->Fill();
			_hTree->Fill();
			_hTaggedTree->Fill();
			ClearVariables();
			vector< Vertex * > unknown = reco.GetUnknownVertexes();
			for (unsigned int i = 0; i < unknown.size(); i++) 
			{
				Write(unknown[i],i);
			}
			std::cout << "We have " << unknown.size() << " unknown vertices\n";
			_numberOfUnknown = unknown.size();
			_hUntaggedTree->Fill();
			ClearVariables();
			delete tagged;
	
		}
		catch( DataNotAvailableException &e)
		{
			streamlog_out(DEBUG) << "No collection!" << std::endl ;
		}
	}
	void  TrashRecoProcessor::Write(LCEvent * evt, std::vector< Particle > * missed,  vector< ReconstructedParticle * > * recomissed)
	{
		IMPL::LCCollectionVec * reco = new IMPL::LCCollectionVec ( LCIO::RECONSTRUCTEDPARTICLE);
		
		_numberOfMissed = missed->size();
		for (unsigned int i = 0; i < missed->size(); i++) 
		{
			if (recomissed) 
			{
				reco->addElement(new ReconstructedParticleImpl((ReconstructedParticleImpl&)(*recomissed->at(i))));

			}
			_momentumMissed[i] = MathOperator::getModule(missed->at(i).GetMomentum());
			_offsetMissed[i] = missed->at(i).GetOffset();
			_thetaMissed[i] = missed->at(i).GetTheta();
		}
		evt->addCollection(reco, _colMissName);
	}
	void TrashRecoProcessor::WriteTagged(vector< VertexTag * > * tagged)
	{
		vector< ReconstructedParticle * > btracks;
		vector< ReconstructedParticle * > bbartracks;
		_bbarnumber1 = 0;
		_bnumber1 = 0;
		_bcharge = 0;
		_bbarcharge = 0;
		int nbvtx = 0;
		int nbbarvtx = 0;
		for (int i = 0; i < _numberOfTagged; i++) 
		{
			Vertex * vertex = tagged->at(i)->GetVertex();
			ReconstructedParticle * particle = vertex->getAssociatedParticle();
			std::cout << "\tPDG: " << tagged->at(i)->GetInitialPDG() <<'\n';
			PrintParticle(particle);
			for (unsigned int j = 0; j < particle->getParticles().size(); j++) 
			{
				PrintParticle(particle->getParticles()[j]);
				if (tagged->at(i)->GetInitialPDG() > 0) 
				{
					btracks.push_back(particle->getParticles()[j]);
				}
				else 
				{
					bbartracks.push_back(particle->getParticles()[j]);
				}
			}
			Write(tagged->at(i), i);
			float distance = MathOperator::getDistance(tagged->at(i)->GetVertex()->getPosition(), _primary->getPosition());
			int particleNumber  = particle->getParticles().size();
			if (tagged->at(i)->GetInitialPDG()  == -5) 
			{
				nbbarvtx++;
				_bbarchimean += tagged->at(i)->GetVertex()->getChi2();
				_bbarprobmean += tagged->at(i)->GetVertex()->getProbability();
				_bbarmomentum += MathOperator::getModule(particle->getMomentum());
				_bbarnumber1 += particleNumber;
				_bbarcharge += particle->getCharge();
				if (_bbarIPdistance < 0.0) 
				{
					_bbarIPdistance = distance;
				}
				else 
				{
					_bbarIPdistance = (_bbarIPdistance > distance)? distance : _bbarIPdistance;
				}
			}
			if (tagged->at(i)->GetInitialPDG() == 5) 
			{
				nbvtx++;
				_bchimean += tagged->at(i)->GetVertex()->getChi2();
				_bprobmean += tagged->at(i)->GetVertex()->getProbability();
				_bmomentum += MathOperator::getModule(particle->getMomentum());
				_bnumber1 += particleNumber;
				_bcharge += particle->getCharge();
				if (_bIPdistance < 0.0) 
				{
					_bIPdistance = distance;
				}
				else 
				{
					_bIPdistance = (_bIPdistance > distance)? distance : _bIPdistance;
				}
			}
		}
		_bbarchimean = (nbbarvtx > 0)? _bbarchimean / (float) nbbarvtx : -1.0;
		_bbarprobmean = (nbbarvtx > 0)? _bbarprobmean / (float) nbbarvtx : -1.0;
		_bchimean = (nbvtx > 0)? _bchimean / (float) nbvtx : -1.0;
		_bprobmean = (nbvtx > 0)? _bprobmean / (float) nbvtx : -1.0;

		_bexists = (_bnumber1 > 0)?1:0;
		_bbarexists = (_bbarnumber1 > 0)?1:0;
		_bptmiss = getMissingPt(btracks, tagged, 5);
		_bbarptmiss = getMissingPt(bbartracks, tagged, -5);
		GetBAngles(tagged);
		std::cout << "Total b multiplicity: " << _bnumber1 << ", total bbar multiplicity: " << _bbarnumber1 << '\n';
		std::cout << "Missing pt b : " << _bptmiss << ", Missing pt bbar: " << _bbarptmiss << '\n';

	}

	void TrashRecoProcessor::WriteTaggedCollection(LCEvent * evt, vector< VertexTag * > * tags)
	{
		IMPL::LCCollectionVec * reco = new IMPL::LCCollectionVec ( LCIO::VERTEX);
		for (int i = 0; i < tags->size(); i++) 
		{
			Vertex * vertex = tags->at(i)->GetVertex();
			reco->addElement(new VertexImpl((VertexImpl&)(*vertex)));

		}
		evt->addCollection(reco, _colTagName);
	}

	void TrashRecoProcessor::GetBAngles(vector< VertexTag * > * tags)
	{
		for (unsigned int i = 0; i < tags->size(); i++) 
		{
			 VertexTag * tag = tags->at(i);
			 const double * vertex = MathOperator::toDoubleArray(tag->GetVertex()->getPosition(),3); 
			 vector< float > direction = MathOperator::getDirection(vertex);
			 if (tag->GetInitialPDG() > 0 && _bteta  < 0.0) 
			 {
			 	_bteta = MathOperator::getAngles(direction)[1];
			 }
			 if (tag->GetInitialPDG() < 0 && _bbarteta  < 0.0) 
			 {
			 	_bbarteta = MathOperator::getAngles(direction)[1];
			 }
		}
	}
	void TrashRecoProcessor::Write(vector< Jet * > * jets)
	{
		_numberOfJets = jets->size();
		for (int i = 0; i < _numberOfJets; i++) 
		{
			_btags[i] = jets->at(i)->GetBTag();
			_ctags[i] = jets->at(i)->GetCTag();
			_mcpdg[i] = jets->at(i)->GetMCPDG();
			_nvertices[i] = jets->at(i)->GetNumberOfVertices();
		}
	}
	float TrashRecoProcessor::getMissingPt(vector< ReconstructedParticle * > & bdaugthers, vector< VertexTag * > * tags, int pdg)
	{
		if (bdaugthers.size() < 1 || !tags || tags->size() < 1) 
		{
			return -1.0;
		}
		VertexTag * bvertex = getBvertex(tags, pdg);
		const float * position = bvertex->GetVertex()->getPosition();
		vector< const double * > vectors;
		for (unsigned int i = 0; i < bdaugthers.size(); i++) 
		{
			vectors.push_back(bdaugthers[i]->getMomentum());
		}
		float missing =  MathOperator::getMissingPt(vectors, position);
		return missing;
	}
	VertexTag * TrashRecoProcessor::getBvertex(vector< VertexTag * > * tags, int pdg)
	{
		if (!tags || tags->size() < 1) 
		{
			return NULL;
		}
		float min = 1000.0;
		int number = -1;
		for (unsigned int i = 0; i < tags->size(); i++) 
		{
			double * pos = MathOperator::toDoubleArray(tags->at(i)->GetVertex()->getPosition() , 3);
			float module = MathOperator::getModule(pos);
			delete pos;
			if (module < min && tags->at(i)->GetInitialPDG() == pdg ) 
			{
				number = i;
				min = module;
			}
		}
		return tags->at(number);
	}
	void TrashRecoProcessor::Write (Vertex * vertex, int number)
	{
		const float * position = vertex->getPosition();
		for (int i = 0; i < 3; i++) 
		{
			_coordinates[number][i] = position[i];
		}
		_probability[number] = vertex->getProbability();
		_chi2[number] = vertex->getChi2();
		ReconstructedParticle * particle = vertex->getAssociatedParticle();
		_charge[number] = particle->getCharge();
		_numberOfParticles[number] = particle->getParticles().size();
		for (int j = 0; j < _numberOfParticles[number]; j++)
		{
			ReconstructedParticle * component = particle->getParticles()[j];
			_momentumOfParticles[number][j] = MathOperator::getModule( component->getMomentum());
			_energyOfParticles[number][j] = component->getEnergy();
			_massOfParticles[number][j] = component->getMass();
		}
		float distance = MathOperator::getDistance(_primary->getPosition(), vertex->getPosition() );
		_distanceFromIP[number] = distance;
		std::cout << "d(V_p,V_i) = " << distance << "; P(V_i) = " << vertex->getProbability() << "; Chi^2(V_i) = " << vertex->getChi2() << '\n';
		
	}
	void TrashRecoProcessor::Write (VertexTag * tag, int number)
	{
		if (!tag) 
		{
			return;
		}
		_numberOfTernary += (tag->GetGeneration() == 3)? 1 : 0;
		_numberOfSecondary += (tag->GetGeneration() == 2)? 1 : 0;
		_PDG[number] = tag->GetInitialPDG();
		_angle[number] = tag->GetTruthAngle();
		_generation[number] = tag->GetGeneration();
		_mindistance[number] = tag->GetMinimalDistance();
		Vertex * vertex = tag->GetVertex();
		Write (vertex, number);
	}
	void TrashRecoProcessor::ClearVariables()
	{
		_bbarchimean = 0.0;
		_bchimean = 0.0;
		_bprobmean = 0.0;
		_bbarprobmean = 0.0;
		_bmomentum = 0.0;
		_bbarmomentum = 0.0;
		_bbarteta = -1.0;
		_bteta = -1.0;
		_numberOfMissed = 0;
		_numberOfJets = 0;
		_numberOfTernary = 0;
		_numberOfSecondary = 0;
		_numberOfTagged = 0;
		_numberOfTotal = 0;
		_numberOfUnknown = 0;
		_bbarnumber1 = -1;
		_bnumber1 = -1;
		_bbarexists = 0;
		_bexists = 0;
		_bbarptmiss = -1.0;
		_bptmiss = -1.0;
		_bcharge = -5;
		_bbarcharge = -5;
		_bbarIPdistance = -1.0;
		_bIPdistance = -1.0;
		for (int i = 0; i < MAXV; i++) 
		{
			_offsetMissed[i] = -1.0;
			_thetaMissed[i] = -1.0;
			_momentumMissed[i] = -1.0;
			_btags[i] = -1.0;
			_ctags[i] = -1.0;
			_nvertices[i] = -1;
			_probability[i] = -1.0;
			_numberOfParticles[i] = -1;
			_PDG[i] = 0;
			_angle[i] = 0.0;
			_generation[i] = 0;
			_chi2[i] = -1.0;
			_distanceFromIP[i] = -1.0;
			for (int j = 0; j < MAXV; j++) 
			{
				_energyOfParticles[i][j] = -1.0;
				_momentumOfParticles[i][j] = -1.0;
				_massOfParticles[i][j] = -1.0;
			}
		}
	}
	
	void TrashRecoProcessor::check( LCEvent * evt ) 
	{ 
		// nothing to check here - could be used to fill checkplots in reconstruction processor
		
	}
	
	
	void TrashRecoProcessor::end()
	{ 
		_hfile->cd();
		_hfile->Write();
		_hfile->Close();
	    //   std::cout << "TrashRecoProcessor::end()  " << name() 
	    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
	    // 	    << std::endl ;
	
	}
} /* TTbarAnalysis */

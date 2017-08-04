
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.io.FileWriter;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.FileNotFoundException;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.random.SynchronizedRandomGenerator;

class queueNetwork {

	/**
	Instance Variables
	*/

	public ArrayList<Double> 				arrivals 	= new ArrayList<Double>();
	public ArrayList<Double> 				arrivalsDom	= new ArrayList<Double>();
	public Double 							arrivalsDomSum = 0.0;
	public ArrayList<Double[]> 				serverRates	= new ArrayList<Double[]>();
	public ArrayList<ArrayList<Double[]>> 	transitions	= new ArrayList<ArrayList<Double[]>>();
	public Double 							domRate		= 0.0; 
	public Double 							scaleDomRate = 2.0; // Calibrate dominating rate multiplier in uniformization
	public int 								jobCount	= 0;
	public Well1024a						rand 		= new Well1024a(2365862);
	public SynchronizedRandomGenerator		rand2 		= new SynchronizedRandomGenerator(rand);	
	public boolean 							fixedTrans 	= true; // Choose whether network topology is fixed or to be inferred
	public Double 							closeProb2	= 0.7; // Calibrate probability of Campling nodes
	public Double 							missingData = 0.5; // Choose percentage of data to be "missing" on average


	/**
	Methods
	*/

	// Method to form a new chain in the network from prior chain and fixed parameters
	public ArrayList<conf> newChain(ArrayList<conf> oldChain){

		double closeProb = closeProb2; 

		// Build a frame with virtual jumps and auxiliaries
		frame newFrame	= addVirtualJumps(oldChain, closeProb);

		// Clamped forward filtering
		ArrayList<ArrayList<confProb>> forwFilt = buildFilter(newFrame.nodes,newFrame.aux);
		
		// When choosing backwards renormalize the probability vector 
		ArrayList<conf> toReturnNew = new ArrayList<conf>(forwFilt.size());

		// Last state
		toReturnNew.add(forwFilt.get(forwFilt.size()-1).get(0).state);
		// Store parent layers in dummy vectors
		ArrayList<Integer> parents = forwFilt.get(forwFilt.size()-1).get(0).upConnect;
		ArrayList<Double> parentsWeights = forwFilt.get(forwFilt.size()-1).get(0).upConnectProbs;
		int chosenParent = 0;
		double normConstant = 0;
		double p 	= 0;
		double sum	= 0;

		for (int i=forwFilt.size()-2; i>=0; i--){ 

			// Use parent layer to choose next entry in sampled sequence
			for (int j=0; j<parents.size(); j++){
				normConstant += forwFilt.get(i).get(parents.get(j)).stateProb * parentsWeights.get(j);
			}

			p 	= rand2.nextDouble();
			sum	= 0;
			outerloop3:
			for (int j=0; j<parents.size(); j++) {
				sum += forwFilt.get(i).get(parents.get(j)).stateProb * parentsWeights.get(j) / normConstant;
				if (p<sum) {
					chosenParent = parents.get(j);
					break outerloop3;
				}
			}
			toReturnNew.add(forwFilt.get(i).get( chosenParent ).state);
			
			// Update parents
			parents = forwFilt.get(i).get( chosenParent ).upConnect;
			parentsWeights = forwFilt.get(i).get( chosenParent ).upConnectProbs;
			normConstant = 0;
			
		}

		// Reorder return vector backwards
		Collections.reverse(toReturnNew);

		// Remove virtual transitions
		int i = 1;
		while (i<toReturnNew.size()){ 
			if (toReturnNew.get(i).configuration.equals(toReturnNew.get(i-1).configuration)) { 
				toReturnNew.remove(i);
			} else i++;
		}

		return toReturnNew;

	}

	// Method to build the forward filter before backward sampling
	public ArrayList<ArrayList<confProb>> buildFilter(ArrayList<conf> frame, ArrayList<int[]> aux){

		ArrayList<ArrayList<confProb>> toReturn = new ArrayList<ArrayList<confProb>>();
		toReturn.add( new ArrayList<confProb>() );
		
		// Starting Configuration always empty
		toReturn.get(0).add( new confProb()); toReturn.get(0).get(0).state = new conf(); toReturn.get(0).get(0).stateProb = 1;

		for (int i=1; i<frame.size(); i++) { // Skip first one, empty and already added above

			if (frame.get(i).in != null) { // Arrival

				toReturn.add( arrivalInFilter2( toReturn.get(toReturn.size()-1), frame, aux, i) );

			} else if (frame.get(i).out != null) { // Departure

				toReturn.add( departureInFilter2( toReturn.get(toReturn.size()-1), frame, aux, i) );

			} else { // Inner Transition / Virtual
				
				toReturn.add( middleJumpInFilter2( toReturn.get(toReturn.size()-1), frame, aux, i) );

			}

		}

		return(toReturn);

	}

	// Method for arrivals BP
	public ArrayList<confProb> arrivalInFilter2(ArrayList<confProb> currentConfs, ArrayList<conf> frame, ArrayList<int[]> aux, int i){

		ArrayList<confProb> toReturn = new ArrayList<confProb>();
		job dummyJob 		= new job( frame.get(i).in.identifier, frame.get(i).in.type );
		confProb dummyConfProb;

		if(currentConfs.size()==0) {
			return(toReturn);
		}

		// If we have a close node we can propagate beliefs quick and easy!
		if (aux.get(i-1)[0] != -2){
			
			for(confProb current : currentConfs){
				dummyConfProb = new confProb();
				dummyConfProb.stateProb 	= current.stateProb; // does not change, arrival likelihood equal to all configurations
				dummyConfProb.upConnect.add(currentConfs.indexOf(current)); // Keep track of origin
				dummyConfProb.state 		= new conf();
				dummyConfProb.state.time 	= frame.get(i).time;
				dummyConfProb.state.in 		= dummyJob;
				// Copy previousConf configurations
				for (int j=0; j<serverRates.size(); j++) {
					dummyConfProb.state.configuration.get(j).addAll(current.state.configuration.get(j));
				}
				// Add the new job
				dummyConfProb.state.configuration.get(aux.get(i-1)[2]).add(dummyJob);
				dummyConfProb.upConnectProbs.add(transitions.get(dummyJob.type).get(0)[aux.get(i-1)[2]]*arrivals.get(dummyJob.type)/domRate); 	// Keep track of origin Prob

				// It's certainly unique
				toReturn.add(dummyConfProb);
			}

			return(toReturn);

		} 

		// Open node...
		for(confProb current : currentConfs){
			int index = 0;
			for (double p : transitions.get(dummyJob.type).get(0)) {
				if (p>0){
					dummyConfProb 				= new confProb();
					dummyConfProb.stateProb 	= current.stateProb * p; 
					dummyConfProb.upConnect.add(currentConfs.indexOf(current)); // Keep track of origin
					dummyConfProb.state 		= new conf();
					dummyConfProb.state.time 	= frame.get(i).time;
					dummyConfProb.state.in 		= dummyJob;
					// Copy previousConf configurations
					for (int j=0; j<serverRates.size(); j++) {
						dummyConfProb.state.configuration.get(j).addAll(current.state.configuration.get(j));
					}
					// Add the new job
					dummyConfProb.state.configuration.get(index).add(dummyJob);
					dummyConfProb.upConnectProbs.add(p * arrivals.get(dummyJob.type) / domRate); 	// Keep track of origin Prob
					// It's certainly unique
					toReturn.add(dummyConfProb);
				}
				index++;
			}
		}

		// No need to normalise
		return(toReturn);

	}

	// Method for departure jumps in chain 
	public ArrayList<confProb> departureInFilter2(ArrayList<confProb> currentConfs, ArrayList<conf> frame, ArrayList<int[]> aux, int i){

		ArrayList<confProb> toReturn = new ArrayList<confProb>();
		confProb dummyConfProb;
		confProb previousConf;

		if(currentConfs.size()==0) {
			return(toReturn);
		}

		// Get state probability vector 
		Double[] probsVectorState = new Double[currentConfs.size()];
		for(int j=0; j<probsVectorState.length; j++) probsVectorState[j] = currentConfs.get(j).stateProb;

		// What servers I need to look at? 
		ArrayList<Integer> exitNodes = new ArrayList<Integer>();
		for(int k=0; k<serverRates.size(); k++){
			if (transitions.get(frame.get(i).out.type).get(k+1)[serverRates.size()]>0.0 ) exitNodes.add(k);
		}
		// Correct if partial observation exists!
		if (frame.get(i).server > 0) {
			exitNodes.clear();
			exitNodes.add(frame.get(i).server-1);
		}
		double[] exitNodesProbs 	= new double[exitNodes.size()];
		for(int k=0; k<exitNodes.size(); k++){
			exitNodesProbs[k] = transitions.get(frame.get(i).out.type).get(exitNodes.get(k)+1)[serverRates.size()] * serverRates.get(exitNodes.get(k))[frame.get(i).out.type] / domRate;
		}

		// Rescale prob-vector and account for closed auxiliaries
		boolean dummyBol = false;
		double normSum = 0.0;
		for(int st=0; st<probsVectorState.length; st++){
			dummyBol = false;
			for (Integer k : exitNodes){
				if (currentConfs.get(st).state.configuration.get(k).size() > 0 && 
					currentConfs.get(st).state.configuration.get(k).get(0).identifier == frame.get(i).out.identifier &&
					(aux.get(i-1)[0]==-2 || aux.get(i-1)[1] == k)){
					dummyBol = true;
				}
			}
			if(!dummyBol){
				probsVectorState[st] = 0.0;
			}
			normSum += probsVectorState[st];
		}

		if(normSum==0.0) {
			return(toReturn);
		}

		for(int j=0; j<probsVectorState.length; j++) probsVectorState[j] = probsVectorState[j]/normSum;

		// If closed auxiliary...
		if (aux.get(i-1)[0] != -2){

			// New configurations from probability vector (whenever not 0 likelihood)
			for(int j=0; j<probsVectorState.length; j++){
				if (probsVectorState[j]>0.0){
					dummyConfProb = new confProb();
					previousConf = currentConfs.get(j);
					dummyConfProb.stateProb = probsVectorState[j]; // remains the same
					dummyConfProb.upConnect.add(j); // Keep track of origin
					dummyConfProb.upConnectProbs.add( exitNodesProbs[ exitNodes.indexOf(aux.get(i-1)[1]) ]); 
					dummyConfProb.state 		= new conf();
					dummyConfProb.state.time 	= frame.get(i).time;
					dummyConfProb.state.out 	= previousConf.state.configuration.get(aux.get(i-1)[1]).get(0);
					dummyConfProb.state.server 	= frame.get(i).server;
					// Copy previousConf configurations
					for (int k=0; k<serverRates.size(); k++) {
						dummyConfProb.state.configuration.get(k).addAll(previousConf.state.configuration.get(k));
					}
					// Remove the job
					dummyConfProb.state.configuration.get(aux.get(i-1)[1]).remove(0);
			 		
			 		// It has to be unique
			 		toReturn.add(dummyConfProb);
				}

			}

			return(toReturn);

		}

		int depQueue = 0;
		double normProb = 0.0;
		for(int j=0; j<probsVectorState.length; j++){
			if (probsVectorState[j]>0.0){
				dummyConfProb = new confProb();
				previousConf = currentConfs.get(j);
				for(int k=0; k<exitNodes.size(); k++){
					if(	previousConf.state.configuration.get(exitNodes.get(k)).size() > 0 && 
						previousConf.state.configuration.get(exitNodes.get(k)).get(0).identifier == frame.get(i).out.identifier){
						depQueue = exitNodes.get(k);
						dummyConfProb.stateProb = probsVectorState[j] * exitNodesProbs[k];
						normProb += dummyConfProb.stateProb;
						dummyConfProb.upConnect.add(j); 
						dummyConfProb.upConnectProbs.add(exitNodesProbs[k]); 
					}
				}
				dummyConfProb.state 		= new conf();
				dummyConfProb.state.time 	= frame.get(i).time;
				dummyConfProb.state.out 	= previousConf.state.configuration.get(depQueue).get(0);
				dummyConfProb.state.server 	= frame.get(i).server;
				// Copy previousConf configurations
				for (int k=0; k<serverRates.size(); k++) {
					dummyConfProb.state.configuration.get(k).addAll(previousConf.state.configuration.get(k));
				}
				// Remove the job
				dummyConfProb.state.configuration.get(depQueue).remove(0);
			 		
		 		//Check if already there before adding any state to Filter
				boolean notThere = true; 
				outerloop:
				for (int k=0; k<toReturn.size(); k++) {
					if ( dummyConfProb.state.configuration.equals(toReturn.get(k).state.configuration ) ){
						toReturn.get(k).stateProb += dummyConfProb.stateProb;
						if(	toReturn.get(k).upConnect.indexOf(dummyConfProb.upConnect.get(0)) == -1 ){
							toReturn.get(k).upConnect.add(dummyConfProb.upConnect.get(0));
							toReturn.get(k).upConnectProbs.add(dummyConfProb.upConnectProbs.get(0));
						}
						notThere = false;
						break outerloop;
					}
				}
				if (notThere) toReturn.add(dummyConfProb);
			}

		}

		//Normalise and return
		for (int j=0; j<toReturn.size(); j++) {
			toReturn.get(j).stateProb = toReturn.get(j).stateProb/normProb;
		}

		return(toReturn);

	}

	// Method for inner jumps in chain 
	public ArrayList<confProb> middleJumpInFilter2(ArrayList<confProb> currentConfs, ArrayList<conf> frame, ArrayList<int[]> aux, int i){

		ArrayList<confProb> toReturn = new ArrayList<confProb>();
		confProb dummyConfProb;
		confProb previousConf;

		if(currentConfs.size()==0) {
			return(toReturn);
		}

		// Get state probability vector
		Double[] probsVectorState = new Double[currentConfs.size()];
		for(int j=0; j<probsVectorState.length; j++) probsVectorState[j] = currentConfs.get(j).stateProb;

		// For later...
		double normConst 	= 0.0;
		int dummyInt 		= 0;
		double weight 		= 0.0;

		// Treat auxiliary cases separately, 
		if (aux.get(i-1)[0] >=0 ){ // So we have a closed node, with a transition (not virtual). Also, served and server redundant

			for(int j=0; j<probsVectorState.length; j++) {
				if ( currentConfs.get(j).state.configuration.get(aux.get(i-1)[1]).size() == 0 ||
					aux.get(i-1)[0] != currentConfs.get(j).state.configuration.get(aux.get(i-1)[1]).get(0).identifier ) {
					probsVectorState[j] = 0.0;
				}
				normConst += probsVectorState[j];
			}

			if(normConst==0.0) {
				return(toReturn);
			}
			for(int j=0; j<probsVectorState.length; j++) probsVectorState[j] = probsVectorState[j] / normConst;

			// Propagate belief...
			for(int j=0; j<probsVectorState.length; j++){
				if (probsVectorState[j]>0.0){
					dummyConfProb 	= new confProb();
					previousConf 	= currentConfs.get(j);
					dummyConfProb.stateProb = probsVectorState[j]; 	// all unique		
					dummyConfProb.upConnect.add(j); 	// Keep track of origin
					dummyConfProb.state 		= new conf();
					dummyConfProb.state.time 	= frame.get(i).time;
					dummyConfProb.state.server 	= frame.get(i).server;
					dummyConfProb.state.served 	= frame.get(i).served;
					// Copy previousConf configurations
					for (int k=0; k<serverRates.size(); k++) {
						dummyConfProb.state.configuration.get(k).addAll(previousConf.state.configuration.get(k));
					}
					// Keep track of parent layer prob
					dummyConfProb.upConnectProbs.add( transitions.get(previousConf.state.configuration.get(aux.get(i-1)[1]).get(0).type).get(aux.get(i-1)[1]+1)[aux.get(i-1)[2]] * 
						serverRates.get(aux.get(i-1)[1])[previousConf.state.configuration.get(aux.get(i-1)[1]).get(0).type] / domRate );	
					// Add and remove
					dummyConfProb.state.configuration.get(aux.get(i-1)[2]).add(previousConf.state.configuration.get(aux.get(i-1)[1]).get(0));
					dummyConfProb.state.configuration.get(aux.get(i-1)[1]).remove(0);
					
					// All unique with closed transition node
				 	toReturn.add(dummyConfProb);

				}
			}	

			return(toReturn);

		} 

		if (aux.get(i-1)[0]==-1){ // Now we have a closed virtual jump, but could be ANY VIRTUAL... served and server cannot exist

			// Propagate belief...
			for(int j=0; j<probsVectorState.length; j++){

				dummyConfProb 	= new confProb();
				previousConf 	= currentConfs.get(j);
				// Weight
				weight = (1-arrivalsDomSum);
				for (int queueIdx = 1; queueIdx <=serverRates.size(); queueIdx++){
					if (previousConf.state.configuration.get(queueIdx-1).size() > 0 ) {
						weight -= serverRates.get(queueIdx-1)[ previousConf.state.configuration.get(queueIdx-1).get(0).type] / domRate;
					}
				}
				dummyConfProb.stateProb = probsVectorState[j] * weight; 	// times virtual likelihood		
				dummyConfProb.upConnect.add(j); 	// Keep track of origin
				normConst		+= dummyConfProb.stateProb;
				dummyConfProb.state 		= new conf();
				dummyConfProb.state.time 	= frame.get(i).time;
				// Copy previousConf configurations
				for (int k=0; k<serverRates.size(); k++) {
					dummyConfProb.state.configuration.get(k).addAll(previousConf.state.configuration.get(k));
				}
				dummyConfProb.upConnectProbs.add(weight);

				// It's virtual... so all unique
				toReturn.add(dummyConfProb);

			}

		} else if (frame.get(i).served > 0){ //account for partial observation if it exist... 

			for ( confProb current : currentConfs ) {
				int queueIdx = frame.get(i).server; 
				if (current.state.configuration.get(queueIdx-1).size() > 0 && 
					current.state.configuration.get(queueIdx-1).get(0).identifier == frame.get(i).served) {
					// Loop through destination queues
					for (int queueIdx2 = 1; queueIdx2 <=serverRates.size(); queueIdx2++) {
						double p = transitions.get( current.state.configuration.get(queueIdx-1).get(0).type ).get(queueIdx)[queueIdx2-1];
						if(p>0){ 
							dummyConfProb 				= new confProb();
							dummyConfProb.stateProb 	= p * current.stateProb *  serverRates.get(queueIdx-1)[ current.state.configuration.get(queueIdx - 1).get(0).type ]  / domRate;
							dummyConfProb.upConnect.add(currentConfs.indexOf(current)); // Keep track of origin
							dummyConfProb.upConnectProbs.add(p *  serverRates.get(queueIdx-1)[ current.state.configuration.get(queueIdx - 1).get(0).type ]  / domRate);
							normConst 					+= dummyConfProb.stateProb;
							dummyConfProb.state 		= new conf();
							dummyConfProb.state.time 	= frame.get(i).time;
							dummyConfProb.state.server 	= frame.get(i).server;
							dummyConfProb.state.served 	= frame.get(i).served;
							// Copy current configurations
							for (int j=0; j<serverRates.size(); j++) {
								dummyConfProb.state.configuration.get(j).addAll(current.state.configuration.get(j));
							}
							// Add and remove the job
							dummyConfProb.state.configuration.get(queueIdx2-1).add(current.state.configuration.get(queueIdx - 1).get(0));
							dummyConfProb.state.configuration.get(queueIdx-1).remove(0);
							//Check if already there before adding any state to Filter
							boolean notThere = true; 
							outerloop:
							for (int j=0; j<toReturn.size(); j++) {
								if ( dummyConfProb.state.configuration.equals(toReturn.get(j).state.configuration) ){
									toReturn.get(j).stateProb += dummyConfProb.stateProb;
									if(	toReturn.get(j).upConnect.indexOf(dummyConfProb.upConnect.get(0)) == -1 ){
										toReturn.get(j).upConnect.add(dummyConfProb.upConnect.get(0));
										toReturn.get(j).upConnectProbs.add(dummyConfProb.upConnectProbs.get(0));
									}
									notThere = false;
									break outerloop;
								}
							}
							if (notThere) toReturn.add(dummyConfProb); 
						}
					}
				}
			}

		} else{ // And now... if the node is open!!! 

			double probsForNoJump = 0; 
			for ( confProb current : currentConfs ) {
				probsForNoJump = 1 - arrivalsDomSum;
				for (int queueIdx = 1; queueIdx <=serverRates.size(); queueIdx++){
					if (current.state.configuration.get(queueIdx-1).size() > 0 ) {
						probsForNoJump	-= serverRates.get(queueIdx-1)[ current.state.configuration.get(queueIdx-1).get(0).type ]  / domRate;
						// Loop through destination queues
						for (int queueIdx2 = 1; queueIdx2 <=serverRates.size(); queueIdx2++) {
							double p = transitions.get( current.state.configuration.get(queueIdx-1).get(0).type ).get(queueIdx)[queueIdx2-1];
							if(p>0){ 
								dummyConfProb 				= new confProb();
								dummyConfProb.stateProb 	= missingData * p * current.stateProb *  serverRates.get(queueIdx-1)[ current.state.configuration.get(queueIdx - 1).get(0).type ]  / domRate;
								dummyConfProb.upConnect.add(currentConfs.indexOf(current)); // Keep track of origin
								dummyConfProb.upConnectProbs.add(missingData * p *  serverRates.get(queueIdx-1)[ current.state.configuration.get(queueIdx - 1).get(0).type ]  / domRate);
								normConst 					+= dummyConfProb.stateProb;
								dummyConfProb.state 		= new conf();
								dummyConfProb.state.time 	= frame.get(i).time;
								// Copy current configurations
								for (int j=0; j<serverRates.size(); j++) {
									dummyConfProb.state.configuration.get(j).addAll(current.state.configuration.get(j));
								}
								// Add and remove the job
								dummyConfProb.state.configuration.get(queueIdx2-1).add(current.state.configuration.get(queueIdx - 1).get(0));
								dummyConfProb.state.configuration.get(queueIdx-1).remove(0);
								//Check if already there before adding any state to Filter
								boolean notThere = true; 
								outerloop:
								for (int j=0; j<toReturn.size(); j++) {
									if ( dummyConfProb.state.configuration.equals(toReturn.get(j).state.configuration) ){
										toReturn.get(j).stateProb += dummyConfProb.stateProb;
										if(	toReturn.get(j).upConnect.indexOf(dummyConfProb.upConnect.get(0)) == -1 ){
											toReturn.get(j).upConnect.add(dummyConfProb.upConnect.get(0));
											toReturn.get(j).upConnectProbs.add(dummyConfProb.upConnectProbs.get(0));
										}
										notThere = false;
										break outerloop;
									}
								}
								if (notThere) toReturn.add(dummyConfProb); 
							}
						}
					}
				}
				// Virtual event
				dummyConfProb 				= new confProb();
				dummyConfProb.stateProb 	= current.stateProb * probsForNoJump; 					
				dummyConfProb.upConnect.add(currentConfs.indexOf(current)); // Keep track of origin
				dummyConfProb.upConnectProbs.add(probsForNoJump);			 // Keep track of origin Prob
				normConst 					+= dummyConfProb.stateProb;
				dummyConfProb.state 		= new conf();
				dummyConfProb.state.time 	= frame.get(i).time;
				// Copy current configurations
				for (int j=0; j<serverRates.size(); j++) {
					dummyConfProb.state.configuration.get(j).addAll(current.state.configuration.get(j));
				}
				//Check if already there before adding any state
				boolean notThere = true; 
				outerloop:
				for (int j=0; j<toReturn.size(); j++) {
					if ( dummyConfProb.state.configuration.equals(toReturn.get(j).state.configuration) ){
						toReturn.get(j).stateProb += dummyConfProb.stateProb;
						if(	toReturn.get(j).upConnect.indexOf(dummyConfProb.upConnect.get(0)) == -1 ){
							toReturn.get(j).upConnect.add(dummyConfProb.upConnect.get(0));
							toReturn.get(j).upConnectProbs.add(dummyConfProb.upConnectProbs.get(0));
						}
						notThere = false;
						break outerloop;
					}
				}
				if (notThere) toReturn.add(dummyConfProb); 
			}

		}
		
		//Normalise and return
		for (int j=0; j<toReturn.size(); j++) {
			toReturn.get(j).stateProb = toReturn.get(j).stateProb/normConst;
		}

		return(toReturn);

	}

	// Method to update parameters given state sequences
	public void updateRates(ArrayList<ArrayList<conf>> stateSequences){

		// To keep record of times and counts
		double totalTime 				= 0.0;
		for (int i =0; i<stateSequences.size(); i++) totalTime += stateSequences.get(i).get(stateSequences.get(i).size()-1).time;
		double[][] timeJobsAtServers 	= new double[serverRates.size()][arrivals.size()];
		double[] countsArr 				= new double[arrivals.size()];
		double[][] countsAtServers 		= new double[serverRates.size()][arrivals.size()];
		double[][][] transCounts		= new double[arrivals.size()][serverRates.size()+1][serverRates.size()+1];


		for (int i=0; i<arrivals.size(); i++) { // Loop through jobs

			double countAr = 0; 											// For counting arrivals of the job
			double[] countAtServer 		= new double[serverRates.size()]; 	// For counting arrivals of the job at each queue
			double[] timeJobAtServer 	= new double[serverRates.size()];	// For counting times at servers for jobs
			double[][] tranCount		= new double[serverRates.size()+1][serverRates.size()+1];  //For counting transitions
			int dummyInteger			= 0;
			boolean dummyBool			= false;
			
			for(int j=0; j<stateSequences.size();j++){
				for(int k=0; k<stateSequences.get(j).size();k++){
					// Task Arrival at Network, and transition counts from arrival :)
					if(stateSequences.get(j).get(k).in != null && stateSequences.get(j).get(k).in.type == i){ 
						countAr++;
						for (int serv=0; serv<serverRates.size(); serv++) {
							dummyInteger = stateSequences.get(j).get(k).configuration.get(serv).size();
							if (dummyInteger>0 && 
								stateSequences.get(j).get(k).configuration.get(serv).get(dummyInteger-1).identifier ==
								stateSequences.get(j).get(k).in.identifier) tranCount[0][serv]++;
						}
					}
					// Task Arrival at each server (NOT QUEUE) and transition counts 
					for (int serv=0; serv<serverRates.size(); serv++) {
						if(stateSequences.get(j).get(k).configuration.get(serv).size() > 0 &&
							stateSequences.get(j).get(k).configuration.get(serv).get(0).type == i){ 
							// Check it's a new job looking at past
							if(stateSequences.get(j).get(k-1).configuration.get(serv).size() == 0 || 
								stateSequences.get(j).get(k).configuration.get(serv).get(0).identifier !=
								stateSequences.get(j).get(k-1).configuration.get(serv).get(0).identifier) {
								countAtServer[serv]++;
								timeJobAtServer[serv] -= stateSequences.get(j).get(k).time;
							}
							// If it's a departure then add the time stayed in the server (NOT QUEUE) and count where transitioned
							if(stateSequences.get(j).get(k+1).configuration.get(serv).size() == 0 || 
								stateSequences.get(j).get(k).configuration.get(serv).get(0).identifier !=
								stateSequences.get(j).get(k+1).configuration.get(serv).get(0).identifier) {
								timeJobAtServer[serv] += stateSequences.get(j).get(k+1).time;
								// Where did it go?
								dummyBool = true;
								for (int serv2=0; serv2<serverRates.size(); serv2++) {
									dummyInteger = stateSequences.get(j).get(k+1).configuration.get(serv2).size();
									if (dummyInteger>0 && 
										stateSequences.get(j).get(k).configuration.get(serv).get(0).identifier ==
										stateSequences.get(j).get(k+1).configuration.get(serv2).get(dummyInteger-1).identifier) {
											tranCount[serv+1][serv2]++;
											dummyBool = false;
									}
								}
								// If out of network?
								if (dummyBool) tranCount[serv+1][serverRates.size()]++;
							}
						}
					}
				}
			}
			countsArr[i] = countAr;
			for (int serv=0; serv<serverRates.size(); serv++) {
				countsAtServers[serv][i] = countAtServer[serv];
				timeJobsAtServers[serv][i] = timeJobAtServer[serv];
				for (int serv2=0; serv2<serverRates.size()+1; serv2++) {
					transCounts[i][serv][serv2] = tranCount[serv][serv2];
				}
			}
			for (int serv2=0; serv2<serverRates.size()+1; serv2++) {
				transCounts[i][serverRates.size()][serv2] = tranCount[serverRates.size()][serv2];
			}

		}

		// POSE A CONTRAINED OPTIMIZATION PROBLEM!
		for (int i=0; i<arrivals.size(); i++) {
			GammaDistribution myGamma = new GammaDistribution(rand2,countsArr[i],1/totalTime);
			arrivals.set(i, myGamma.sample()); 

			double truncUnif = 0;
			for (int serv=0; serv<serverRates.size(); serv++) {
				myGamma = new GammaDistribution(rand2,countsAtServers[serv][i],1/timeJobsAtServers[serv][i]);

				if (serv == 0){
					truncUnif = rand2.nextDouble() * myGamma.cumulativeProbability(serverRates.get(serv+1)[i]);
				} else if (serv == serverRates.size() - 1){
					truncUnif = myGamma.cumulativeProbability(serverRates.get(serv-1)[i]) + 
						rand2.nextDouble() * ( 1 - myGamma.cumulativeProbability(serverRates.get(serv-1)[i]) );
				} else {
					truncUnif = myGamma.cumulativeProbability(serverRates.get(serv-1)[i]) + 
						rand2.nextDouble() * ( myGamma.cumulativeProbability(serverRates.get(serv+1)[i]) - myGamma.cumulativeProbability(serverRates.get(serv-1)[i]) );				
				}

				serverRates.get(serv)[i] = myGamma.inverseCumulativeProbability(truncUnif);

			}
		}

		// ONLY USED IF WANT TO DRAW INFERENCE ON TOPOLOGY
		// Probabilities go Gamma -> Ratios, and that is Dirichlet
		// if (!fixedTrans){
		// 	for (int job=0; job<arrivals.size(); job++) {
		// 		double[] runSum = new double[serverRates.size()+1];
		// 		for (int i=0; i<=serverRates.size() ; i++) {
		// 			for (int j=0; j<=serverRates.size() ; j++) {
		// 				if(transCounts[job][i][j]>0.0){ //If it was 0, the that route is user blocked, no update and remains at 0
		// 					GammaDistribution myGamma = new GammaDistribution(rand2,transCounts[job][i][j],1);
		// 					transitions.get(job).get(i)[j] = myGamma.sample();
		// 					runSum[i] += transitions.get(job).get(i)[j];
		// 				}
		// 			}
		// 			//Ratios
		// 			for (int j=0; j<=serverRates.size() ; j++) {
		// 				transitions.get(job).get(i)[j] = transitions.get(job).get(i)[j]/runSum[i];
		// 			}
		// 		}	
		// 	}
		// }

		// Work out Dominating Rate for Uniformization
		domRate	= 0.0;
		for (Double i : arrivals){
			domRate += i;
		}
		for (Double[] i : serverRates){
			domRate += Collections.max(Arrays.asList(i));
		}
		// Scale for domRate, nees to be strictly greater than sum of max rates!
		domRate *= scaleDomRate;
		// Add dominated arrival rates
		for (int i=0; i<arrivals.size(); i++ ) {
			arrivalsDom.set(i,arrivals.get(i)/domRate);
		}
		System.out.println("\nNetwork Updated, the strictly greater dominating rate is: " + domRate + "\n");
		//Sum them over to avoid cumbersome operations later
		arrivalsDomSum = 0.0;
		for (int j=0; j<arrivalsDom.size(); j++) {
			arrivalsDomSum += arrivalsDom.get(j);
		}

	}

	// Arbitrary Initial sequences in Bottleneck Network
	// Applies random branch reduction to handle massive state space
	public ArrayList<conf> initFFBS(ArrayList<obs> dummyObs){

		// Departure order is helpful for our case... build first network by heuristics
		ArrayList<Integer> depOrder = new ArrayList<Integer>();
		for (int i=0; i<dummyObs.size(); i++){
			if (dummyObs.get(i).out>0) depOrder.add(dummyObs.get(i).out);
		}

		// Create a sequence of families of states in the filter
		ArrayList<ArrayList<confProb>> forwFilt = initFilter(dummyObs,1000,1000,depOrder);
		// Many will be empty, rejection until at least one possible backward sequence
		while ( forwFilt.get(forwFilt.size()-1).size() == 0) {
			forwFilt = initFilter(dummyObs,1000,1000,depOrder);
		}

		// Create return vector, and sequence of families of states in the filter
		ArrayList<conf> toReturn = new ArrayList<conf>(forwFilt.size());
		// Now goes backward sampling; do in wrong order, then reorder
		toReturn.add(forwFilt.get(forwFilt.size()-1).get(0).state);

		// Store parent layer in dummy vector
		ArrayList<Integer> parents = forwFilt.get(forwFilt.size()-1).get(0).upConnect;
		ArrayList<Double> parentsWeights = forwFilt.get(forwFilt.size()-1).get(0).upConnectProbs;
		int chosenParent = 0;
		double normConstant = 0;
		double p 	= 0;
		double sum	= 0;

		for (int i=forwFilt.size()-2; i>=0; i--){ 

			// Use parent layer to randomly choose next entry in sampled sequence
			for (int j=0; j<parents.size(); j++){
				normConstant += forwFilt.get(i).get(parents.get(j)).stateProb * parentsWeights.get(j);
			}
			p 	= rand2.nextDouble();
			sum	= 0;
			outerloop3:
			for (int j=0; j<parents.size(); j++) {
				sum += forwFilt.get(i).get(parents.get(j)).stateProb * parentsWeights.get(j) / normConstant;
				if (p<sum) {
					chosenParent = parents.get(j);
					break outerloop3;
				}
			}
			toReturn.add(forwFilt.get(i).get( chosenParent ).state);

			// Update parents
			parents = forwFilt.get(i).get( chosenParent ).upConnect;
			parentsWeights = forwFilt.get(i).get( chosenParent ).upConnectProbs;
			normConstant = 0;
			
		}

		// Reorder return vector backwards
		Collections.reverse(toReturn);
		return toReturn;

	}

	// Method to build the forward filter before backward sampling within Initializiation
	public ArrayList<ArrayList<confProb>> initFilter(ArrayList<obs> dummyObs, int maxArr, int maxInner, ArrayList<Integer> depOrder){

		ArrayList<ArrayList<confProb>> toReturn = new ArrayList<ArrayList<confProb>>();
		//First make some virtual jumps, and go filtering 
		toReturn.add( new ArrayList<confProb>() );
		// Starting Configuration always empty
		toReturn.get(0).add( new confProb()); toReturn.get(0).get(0).state = new conf(); toReturn.get(0).get(0).stateProb = 1;
		ArrayList<Integer> controlOrders = new ArrayList<Integer>();
		controlOrders.addAll(depOrder);

		for (int i=0; i<dummyObs.size(); i++) {

			// Need 1 virtual jump before departure in BOTTLENECK QUEUE
			if (dummyObs.get(i).in >-1) {

				//Arrival Jump
				toReturn.add( arrivalInFilterInit( toReturn.get(toReturn.size()-1), dummyObs, i , maxArr, depOrder) );

			} else if (dummyObs.get(i).served >-1) {

				// Retain departure orders with missing info!
				while (controlOrders.get(0) != dummyObs.get(i).served){
					// Middle jump in advance
					toReturn.add( priorInnerInit(toReturn.get(toReturn.size()-1), dummyObs, i, controlOrders.get(0)) );
					controlOrders.remove(0);
				}
				toReturn.add( observedTransitionInit(toReturn.get(toReturn.size()-1), dummyObs, i) );
				controlOrders.remove(0);
				
			} else if (dummyObs.get(i).out >-1) {

				if (controlOrders.size()>0 && controlOrders.get(0) == dummyObs.get(i).out){
					// Middle jump in advance
					toReturn.add( middleJumpInFilterInit(toReturn.get(toReturn.size()-1), dummyObs, i, maxInner) );
					controlOrders.remove(0);
				}
				// Jump to departure
				toReturn.add( departureInFilterInit(toReturn.get(toReturn.size()-1), dummyObs, i) );
				
			}

		}

		return(toReturn);

	}

	// Method for arrivals when initializing chain (random chopping)
	public ArrayList<confProb> arrivalInFilterInit(ArrayList<confProb> currentConfs, ArrayList<obs> dummyObs, int i,
		int maxSize, ArrayList<Integer> depOrder){

		ArrayList<confProb> toReturn = new ArrayList<confProb>();
		job dummyJob = new job( dummyObs.get(i).in, dummyObs.get(i).inType-1 );
		confProb dummyConfProb;
		double probsSum = 0;
		for ( confProb previousConf : currentConfs ) {
			// It's an arrival, probs given by transition matrix
			int index = 0;
			for (double p : transitions.get(dummyJob.type).get(0)) {
				if(p>0 &&
					// Only if no job there that departs after it
					(previousConf.state.configuration.get(index).size() == 0 ||	
						depOrder.indexOf(dummyJob.identifier)> depOrder.indexOf(previousConf.state.configuration.get(index).get(previousConf.state.configuration.get(index).size()-1).identifier)) ){
					dummyConfProb 				= new confProb();
					dummyConfProb.stateProb 	= p * previousConf.stateProb;
					dummyConfProb.upConnect.add(currentConfs.indexOf(previousConf)); // Keep track of origin
					dummyConfProb.upConnectProbs.add(p); // Keep track of origin Prob
					probsSum					+= dummyConfProb.stateProb;
					dummyConfProb.state 		= new conf();
					dummyConfProb.state.time 	= dummyObs.get(i).time;
					dummyConfProb.state.in 		= dummyJob;
					// Copy previousConf configurations
					for (int j=0; j<serverRates.size(); j++) {
						dummyConfProb.state.configuration.get(j).addAll(previousConf.state.configuration.get(j));
					}
					// Add the new job
					dummyConfProb.state.configuration.get(index).add(dummyJob);
					// Add to the filter
					toReturn.add(dummyConfProb); 
				}
				index++;
			} 
		}
		// Take maxSize random states and wish for it to be Ok! Increase values if cannot get backward sampling
		if (toReturn.size()>maxSize){
			toReturn = randChopNew(toReturn,maxSize);
		} else {
			for (int j=0; j<toReturn.size(); j++) {
				toReturn.get(j).stateProb = toReturn.get(j).stateProb/probsSum;
			}
		}

		return(toReturn);

	}

	// Method for inner jumps when initializing chain
	public ArrayList<confProb> observedTransitionInit(ArrayList<confProb> currentConfs, ArrayList<obs> dummyObs, int i){

		ArrayList<confProb> toReturn = new ArrayList<confProb>();
		confProb dummyConfProb;
		double probsSum = 0;
		// Server is fixed...
		int queueIdx = dummyObs.get(i).server;

		for ( confProb previousConf : currentConfs ) {
			// Only if the queue is populated
			if (previousConf.state.configuration.get(queueIdx-1).size() > 0 &&
				// Just look for the job departing, otherwise too many options
				previousConf.state.configuration.get(queueIdx-1).get(0).identifier == dummyObs.get(i).served ) {
				// Loop through destination queues
				for (int queueIdx2 = 1; queueIdx2 <=serverRates.size(); queueIdx2++) {
					double p = transitions.get( previousConf.state.configuration.get(queueIdx-1).get(0).type ).get(queueIdx)[queueIdx2-1];
					if(p>0){ 
						dummyConfProb 				= new confProb();
						dummyConfProb.stateProb 	= p * previousConf.stateProb *  serverRates.get(queueIdx-1)[ previousConf.state.configuration.get(queueIdx - 1).get(0).type ]  / domRate;
						dummyConfProb.upConnect.add(currentConfs.indexOf(previousConf)); // Keep track of origin
						dummyConfProb.upConnectProbs.add(p *  serverRates.get(queueIdx-1)[ previousConf.state.configuration.get(queueIdx - 1).get(0).type ]  / domRate);
						probsSum 					+= dummyConfProb.stateProb;
						dummyConfProb.state 		= new conf();
						dummyConfProb.state.time 	= dummyObs.get(i).time;
						dummyConfProb.state.served 	= dummyObs.get(i).served;
						dummyConfProb.state.server 	= dummyObs.get(i).server;
						// Copy previousConf configurations
						for (int j=0; j<serverRates.size(); j++) {
							dummyConfProb.state.configuration.get(j).addAll(previousConf.state.configuration.get(j));
						}
						// Add and remove the job
						dummyConfProb.state.configuration.get(queueIdx2-1).add(previousConf.state.configuration.get(queueIdx - 1).get(0));
						dummyConfProb.state.configuration.get(queueIdx-1).remove(0);
						// All unique
						toReturn.add(dummyConfProb); 
					}
				}
			}
		}

		for (int j=0; j<toReturn.size(); j++) {
			toReturn.get(j).stateProb = toReturn.get(j).stateProb/probsSum;
		}

		return(toReturn);

	}

	// Method for inner jumps when initializing chain
	public ArrayList<confProb> priorInnerInit(ArrayList<confProb> currentConfs, ArrayList<obs> dummyObs, int i,	int idMove){

		ArrayList<confProb> toReturn = new ArrayList<confProb>();
		confProb dummyConfProb;
		double probsSum = 0;
 
		for ( confProb previousConf : currentConfs ) {
			// It's a forced transition without arrivals or departures
			for (int queueIdx = 1; queueIdx <=serverRates.size(); queueIdx++){
			// Only if the queue is populated
				if (previousConf.state.configuration.get(queueIdx-1).size() > 0 &&
					// Just look for the job departing, otherwise too many options
					previousConf.state.configuration.get(queueIdx-1).get(0).identifier == idMove ) {
					// Loop through destination queues
					for (int queueIdx2 = 1; queueIdx2 <=serverRates.size(); queueIdx2++) {
						double p = transitions.get( previousConf.state.configuration.get(queueIdx-1).get(0).type ).get(queueIdx)[queueIdx2-1];
						if(p>0){ 
							dummyConfProb 				= new confProb();
							dummyConfProb.stateProb 	= p * previousConf.stateProb *  serverRates.get(queueIdx-1)[ previousConf.state.configuration.get(queueIdx - 1).get(0).type ]  / domRate;
							dummyConfProb.upConnect.add(currentConfs.indexOf(previousConf)); // Keep track of origin
							dummyConfProb.upConnectProbs.add(p *  serverRates.get(queueIdx-1)[ previousConf.state.configuration.get(queueIdx - 1).get(0).type ]  / domRate);
							probsSum 					+= dummyConfProb.stateProb;
							dummyConfProb.state 		= new conf();
							dummyConfProb.state.time 	= previousConf.state.time + (dummyObs.get(i).time - previousConf.state.time) / 2.0;
							// Copy previousConf configurations
							for (int j=0; j<serverRates.size(); j++) {
								dummyConfProb.state.configuration.get(j).addAll(previousConf.state.configuration.get(j));
							}
							// Add and remove the job
							dummyConfProb.state.configuration.get(queueIdx2-1).add(previousConf.state.configuration.get(queueIdx - 1).get(0));
							dummyConfProb.state.configuration.get(queueIdx-1).remove(0);
							//Check if already there before adding any state to Filter
							boolean notThere = true; 
							outerloop:
							for (int j=0; j<toReturn.size(); j++) {
								if ( dummyConfProb.state.configuration.equals(toReturn.get(j).state.configuration) ){
									toReturn.get(j).stateProb += dummyConfProb.stateProb;
									if(	toReturn.get(j).upConnect.indexOf(dummyConfProb.upConnect.get(0)) == -1 ){
										toReturn.get(j).upConnect.add(dummyConfProb.upConnect.get(0));
										toReturn.get(j).upConnectProbs.add(dummyConfProb.upConnectProbs.get(0));
									}
									notThere = false;
									break outerloop;
								}
							}
							if (notThere) toReturn.add(dummyConfProb); 
						}
					}
				}
			}
		}

		for (int j=0; j<toReturn.size(); j++) {
			toReturn.get(j).stateProb = toReturn.get(j).stateProb/probsSum;
		}

		return(toReturn);

	}

	// Method for inner jumps when initializing chain 
	public ArrayList<confProb> middleJumpInFilterInit(ArrayList<confProb> currentConfs, ArrayList<obs> dummyObs, int i,
		int maxSize){

		ArrayList<confProb> toReturn = new ArrayList<confProb>();
		confProb dummyConfProb;
		double probsSum = 0;
		//double probsForNoJump = 0; 
		for ( confProb previousConf : currentConfs ) {
			// It's a forced transition without arrivals or departures
			//probsForNoJump = arrivalsDomSum;
			for (int queueIdx = 1; queueIdx <=serverRates.size(); queueIdx++){
			// Only if the queue is populated
				if (previousConf.state.configuration.get(queueIdx-1).size() > 0 &&
					// Just look for the job departing, otherwise too many options
					previousConf.state.configuration.get(queueIdx-1).get(0).identifier == dummyObs.get(i).out ) {
					//probsForNoJump	+= serverRates.get(queueIdx-1)[ previousConf.state.configuration.get(queueIdx-1).get(0).type ]  / domRate;
				// Loop through destination queues
					for (int queueIdx2 = 1; queueIdx2 <=serverRates.size(); queueIdx2++) {
						double p = transitions.get( previousConf.state.configuration.get(queueIdx-1).get(0).type ).get(queueIdx)[queueIdx2-1];
						if(p>0){ //Here goes the damm Beam sampler uniform variable to cut this mess down!
							dummyConfProb 				= new confProb();
							dummyConfProb.stateProb 	= p * previousConf.stateProb *  serverRates.get(queueIdx-1)[ previousConf.state.configuration.get(queueIdx - 1).get(0).type ]  / domRate;
							dummyConfProb.upConnect.add(currentConfs.indexOf(previousConf)); // Keep track of origin
							dummyConfProb.upConnectProbs.add(p *  serverRates.get(queueIdx-1)[ previousConf.state.configuration.get(queueIdx - 1).get(0).type ]  / domRate);
							probsSum 					+= dummyConfProb.stateProb;
							dummyConfProb.state 		= new conf();
							dummyConfProb.state.time 	= previousConf.state.time + (dummyObs.get(i).time - previousConf.state.time) / 2.0;
							// Copy previousConf configurations
							for (int j=0; j<serverRates.size(); j++) {
								dummyConfProb.state.configuration.get(j).addAll(previousConf.state.configuration.get(j));
							}
							// Add and remove the job
							dummyConfProb.state.configuration.get(queueIdx2-1).add(previousConf.state.configuration.get(queueIdx - 1).get(0));
							dummyConfProb.state.configuration.get(queueIdx-1).remove(0);
							//Check if already there before adding any state to Filter
							boolean notThere = true; 
							outerloop:
							for (int j=0; j<toReturn.size(); j++) {
								if ( dummyConfProb.state.configuration.equals(toReturn.get(j).state.configuration) ){
									toReturn.get(j).stateProb += dummyConfProb.stateProb;
									if(	toReturn.get(j).upConnect.indexOf(dummyConfProb.upConnect.get(0)) == -1 ){
										toReturn.get(j).upConnect.add(dummyConfProb.upConnect.get(0));
										toReturn.get(j).upConnectProbs.add(dummyConfProb.upConnectProbs.get(0));
									}
									notThere = false;
									break outerloop;
								}
							}
							if (notThere) toReturn.add(dummyConfProb); 
						}
					}
				}
			}
		}

		// Take maxSize random states and wish for it to be Ok! Increase these values if cannot get backward sampling
		if (toReturn.size()>maxSize){
			toReturn = randChopNew(toReturn,maxSize); // Comes normalised
		} else {
			for (int j=0; j<toReturn.size(); j++) {
				toReturn.get(j).stateProb = toReturn.get(j).stateProb/probsSum;
			}
		}

		return(toReturn);

	}

	// Method for inner jumps when initializing chain (no chopping, already restrictive :) )
	public ArrayList<confProb> departureInFilterInit(ArrayList<confProb> currentConfs, ArrayList<obs> dummyObs, int i){

		ArrayList<confProb> toReturn = new ArrayList<confProb>();
		confProb dummyConfProb;
		double probsSum = 0;
		for ( confProb previousConf : currentConfs ) {
			// It's a departure, and so we only consider such transitions
			for (int queueIdx = 1; queueIdx <=serverRates.size(); queueIdx++){
				// Only if the queue is populated and can transition towards a departure
				if (previousConf.state.configuration.get(queueIdx-1).size() > 0 &&
					//Need be the very job that we observe departing
					previousConf.state.configuration.get(queueIdx-1).get(0).identifier == dummyObs.get(i).out &&
					// The beam sampler uniform variable goes inside this conditional!
					transitions.get(previousConf.state.configuration.get(queueIdx-1).get(0).type).get(queueIdx)[serverRates.size()] > 0.0) {
					// So we remove job :)
					double p = transitions.get(previousConf.state.configuration.get(queueIdx-1).get(0).type).get(queueIdx)[serverRates.size()];
					dummyConfProb 				= new confProb();
					dummyConfProb.stateProb 	= p * previousConf.stateProb *  serverRates.get(queueIdx-1)[ previousConf.state.configuration.get(queueIdx - 1).get(0).type ]  / domRate;
					dummyConfProb.upConnect.add(currentConfs.indexOf(previousConf)); // Keep track of origin
					dummyConfProb.upConnectProbs.add(p *  serverRates.get(queueIdx-1)[ previousConf.state.configuration.get(queueIdx - 1).get(0).type ]  / domRate);
					probsSum 					+= dummyConfProb.stateProb;
					dummyConfProb.state 		= new conf();
					dummyConfProb.state.time 	= dummyObs.get(i).time;
					dummyConfProb.state.out 	= previousConf.state.configuration.get(queueIdx-1).get(0);
					if(dummyObs.get(i).server > 0)	dummyConfProb.state.server 	= dummyObs.get(i).server;
					// Copy previousConf configurations
					for (int j=0; j<serverRates.size(); j++) {
						dummyConfProb.state.configuration.get(j).addAll(previousConf.state.configuration.get(j));
					}
					// Remove the job
					dummyConfProb.state.configuration.get(queueIdx-1).remove(0);
					//Check if already there before adding any state to Filter
					boolean notThere = true; 
					outerloop:
					for (int j=0; j<toReturn.size(); j++) {
						if ( dummyConfProb.state.configuration.equals(toReturn.get(j).state.configuration ) ){
							toReturn.get(j).stateProb += dummyConfProb.stateProb;
							if(	toReturn.get(j).upConnect.indexOf(dummyConfProb.upConnect.get(0)) == -1 ){
								toReturn.get(j).upConnect.add(dummyConfProb.upConnect.get(0));
								toReturn.get(j).upConnectProbs.add(dummyConfProb.upConnectProbs.get(0));
							}
							notThere = false;
							break outerloop;
						}
					}
					if (notThere) toReturn.add(dummyConfProb); 
				}
			}
		}
		//Normalise
		for (int j=0; j<toReturn.size(); j++) {
			toReturn.get(j).stateProb = toReturn.get(j).stateProb/probsSum;
		}

		return(toReturn);

	}

	// Method for likelihood based (then randomly) chop in a filter at some point
	public ArrayList<confProb> randChopNew(ArrayList<confProb> currentConfs, int amount){

		ArrayList<confProb> toReturn = new ArrayList<confProb>();
		double[] probsStates = new double[currentConfs.size()];
		for (int j=0; j<currentConfs.size(); j++) {
			probsStates[j] = currentConfs.get(j).stateProb;
		}

		double sum = 0.0;
		double p = 0.0;
		double fullProb = 1.0;

		for (int i = 0; i<amount; i++){
			p = rand2.nextDouble() * fullProb;
			sum	= 0;
			outerloop5:
			for (int j=0; j<probsStates.length; j++) {
				sum += probsStates[j];
				if (p<sum) {
					toReturn.add(currentConfs.get(j));
					fullProb -= probsStates[j];
					probsStates[j] = 0.0;
					break outerloop5;
				}
			}
		}

		// If I do cuts... need to normalise again!
		for (int j=0; j<toReturn.size(); j++) {
			toReturn.get(j).stateProb = toReturn.get(j).stateProb/(1-fullProb);
		}		

		return(toReturn);

	}

	// To read in simulation files and manipulate them for inference purposes
	public ArrayList<obs> readSimulation(String Filename){

		ArrayList<obs> retObs	= new ArrayList<obs>();
		BufferedReader br 		= null;
		String line 			= "";
		String cvsSplitBy 		= ", ";

		try{
			br = new BufferedReader(new FileReader(System.getProperty("user.dir") + "\\Files\\" + Filename));
			int index = 0;
			br.readLine(); // Skip first line
			while ( (line = br.readLine()) != null ) {
				// use comma as separator
				String[] splitLine 	= line.split(cvsSplitBy);
				obs dummyObs		= new obs();
				dummyObs.time 		= Double.parseDouble(splitLine[0]);
				if ( index>0 )		dummyObs.jobsInIDs.addAll(retObs.get(index-1).jobsInIDs);
				if (!splitLine[1].equals("")) {
					dummyObs.in		= Integer.parseInt(splitLine[1]); 
					dummyObs.inType = Integer.parseInt(splitLine[2]);
					dummyObs.jobsInIDs.add(dummyObs.in);
				}
				if (!splitLine[3].equals("")) {
					dummyObs.out 	= Integer.parseInt(splitLine[3]);
					if (rand2.nextDouble() > missingData) dummyObs.server = Integer.parseInt(splitLine[5]);					 
					dummyObs.jobsInIDs.remove( (Integer) dummyObs.out);
				}
				if (!splitLine[4].equals("") && rand2.nextDouble() > missingData) {
					dummyObs.served = Integer.parseInt(splitLine[4]);	
					dummyObs.server = Integer.parseInt(splitLine[5]);	
				}
				retObs.add(dummyObs);
				index++;
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally{
			try {
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		return retObs;

	}

	// To simulate and output a Sequence of Configurations, up to a provided terminal time, And write it!!!
	public void writeSimulationUnbounded(double T, String Filename){
		
		try{
			FileWriter writer = new FileWriter(System.getProperty("user.dir") + Filename);
			writer.append("Time" + ", " + "In" + ", " + "Type" + ", " + "Out" + ", " + "Served" + ", " + "Server" + ", " + "Current Jobs" + "\n");
			conf startConf 	= new conf();
			conf dummyConf;
			startConf.time 	= 0;
			jobCount 		= 0;

			int totJobs = 0;
			int cumJobs = 0;
			while (startConf.time < T || totJobs > 0){

				dummyConf = new conf();
				for (int j=0; j<serverRates.size(); j++) {
					dummyConf.configuration.get(j).addAll(startConf.configuration.get(j));
				}
				update(startConf);

				if(startConf.time >= T && totJobs == 0){
					writer.flush();
					writer.close();
					System.out.println("Simulation Completed");
					if (cumJobs == 0) writeSimulationUnbounded(T,Filename);
					return;					
				}
				if(startConf.in != null){
					writer.append(startConf.time + ", " + startConf.in.identifier+ ", " + (startConf.in.type + 1) + ", " + ", ");
					totJobs = 0;
					cumJobs++;
					for (int i=0; i<serverRates.size(); i++) totJobs += startConf.configuration.get(i).size();
					writer.append(", " + ", " + totJobs + "\n");						
				} else if(startConf.out != null){
					writer.append(startConf.time + ", " + ", " + ", " + startConf.out.identifier + ", ");
					totJobs = 0;
					for (int i=0; i<serverRates.size(); i++) totJobs += startConf.configuration.get(i).size();
					for (int j=0; j<serverRates.size(); j++) {
						if (dummyConf.configuration.get(j).size()>0 && 
							(startConf.configuration.get(j).size() == 0 || 
								!dummyConf.configuration.get(j).get(0).equals(startConf.configuration.get(j).get(0)) )){
							writer.append(", " + (j+1) + ", " + totJobs + "\n");
						}			
					}
				} else {
					for (int j=0; j<serverRates.size(); j++) {
						if (dummyConf.configuration.get(j).size()>0 && 
							(startConf.configuration.get(j).size() == 0 || 
								!dummyConf.configuration.get(j).get(0).equals(startConf.configuration.get(j).get(0)) )){
							writer.append(startConf.time + ", " + ", " + ", " + ", ");
							writer.append(dummyConf.configuration.get(j).get(0).identifier + ", " + (j+1) + ", " + totJobs + "\n");
						}			
					}
				}

			}
			writer.flush();
			writer.close();
			System.out.println("Simulation Completed");
			if (cumJobs == 0) writeSimulationUnbounded(T,Filename);
			return;
		} catch (IOException e){
			e.printStackTrace();
		}

	}

	// To simulate and output a Sequence of Configurations, up to a provided terminal time, AND RETURN IT!
	public ArrayList<conf> doSimulation(double T){
		
		ArrayList<conf> returnSeq = new ArrayList<conf>();
		returnSeq.add(new conf());
		returnSeq.get(0).time 	= 0;

		int totJobs = 0;
		int jobCountDummy = 0;
		outerloop4:
		while (returnSeq.get(returnSeq.size()-1).time < T || totJobs > 0){
			returnSeq.add( updateAlternative( returnSeq.get( returnSeq.size()-1 ) ,  jobCountDummy) );
			if(returnSeq.get(returnSeq.size()-1).time >= T && totJobs == 0){
				returnSeq.remove( returnSeq.size()-1 );
				break outerloop4;					
			}
			if(returnSeq.get(returnSeq.size()-1).in != null){
				totJobs++;				
				jobCountDummy++;
			} else if(returnSeq.get(returnSeq.size()-1).out != null){
				totJobs --;
			}
		}
		// REMOVE VIRTUALS
		int i = 1;
		while (i<returnSeq.size()){ 
			if (returnSeq.get(i).configuration.equals(returnSeq.get(i-1).configuration)) { 
				returnSeq.remove(i);
			} else i++;
		}
		if (jobCountDummy == 0) return (doSimulation(T));
		return returnSeq;

	}

	// Simulate a Realization from the Queue, assume rates and time reference to be in hours
	public void update(conf dummyConf) {
		
		dummyConf.time = dummyConf.time + Math.log(1-rand2.nextDouble())/(-domRate);
		dummyConf.in = null;	
		dummyConf.out = null;	
		/*FCFS Policy*/
		double prob = rand2.nextDouble();
		double sum 	= 0;
		for (int i=0; i<arrivalsDom.size(); i++) {
			sum += arrivalsDom.get(i);
			if (prob < sum) {
				//Transition is an arrival
				jobCount++;
				job newJob = new job(jobCount,i);
				int destServer = probUtils.sampleIndex(transitions.get(i).get(0),rand2);
				dummyConf.configuration.get(destServer).add(newJob);
				dummyConf.in = newJob;
				return;
			}
		}
		for (int i=0; i<serverRates.size(); i++) {
			if (dummyConf.configuration.get(i).size() > 0){
				sum +=  serverRates.get(i)[((job) dummyConf.configuration.get(i).get(0)).type] / domRate;
				if (prob < sum) {
					//Transition from that server
					int destServer = probUtils.sampleIndex(transitions.get(((job) dummyConf.configuration.get(i).get(0)).type).get(i+1),rand2);
					if (destServer == serverRates.size()){
						dummyConf.out = dummyConf.configuration.get(i).get(0);
						dummyConf.configuration.get(i).remove(0);
						return;
					} else {
						dummyConf.configuration.get(destServer).add((job) dummyConf.configuration.get(i).get(0));
						dummyConf.configuration.get(i).remove(0);
						return;
					}
				}
			}
		}

	} 

	// Simulate a Realization from the Queue, assume rates and time reference to be in hours, return new object!
	public conf updateAlternative(conf dummyConf, int jobCountDummy) {
		
		conf returnConf = new conf();
		returnConf.time = dummyConf.time + Math.log(1-rand2.nextDouble())/(-domRate);
		returnConf.in = null;	
		returnConf.out = null;	

		// Copy previousConf configurations
		for (int j=0; j<serverRates.size(); j++) {
			returnConf.configuration.get(j).addAll(dummyConf.configuration.get(j));
		}

		/*FCFS Policy*/

		double prob = rand2.nextDouble();
		double sum 	= 0;
		for (int i=0; i<arrivalsDom.size(); i++) {
			sum += arrivalsDom.get(i);
			if (prob < sum) {
				//Transition is an arrival
				job newJob = new job(jobCountDummy + 1,i);
				int destServer = probUtils.sampleIndex(transitions.get(i).get(0),rand2);
				returnConf.configuration.get(destServer).add(newJob);
				returnConf.in = newJob;
				return returnConf;
			}
		}
		for (int i=0; i<serverRates.size(); i++) {
			if (returnConf.configuration.get(i).size() > 0){
				sum +=  serverRates.get(i)[((job) returnConf.configuration.get(i).get(0)).type] / domRate;
				if (prob < sum) {
					//Transition from that server
					int destServer = probUtils.sampleIndex(transitions.get(((job) returnConf.configuration.get(i).get(0)).type).get(i+1),rand2);
					if (destServer == serverRates.size()){
						returnConf.out = returnConf.configuration.get(i).get(0);
						returnConf.configuration.get(i).remove(0);
						return returnConf;
					} else {
						returnConf.configuration.get(destServer).add((job) returnConf.configuration.get(i).get(0));
						returnConf.configuration.get(i).remove(0);
						return returnConf;
					}
				}
			}
		}

		return returnConf;

	} 

	// Method to construct a frame with virtual jumps taking a current chain of jumps (so ignore observations -> redundant)
	public frame addVirtualJumps(ArrayList<conf> oldChain, double closeProb){

		frame toReturn 	= new frame();
		
		// Add virtual jumps, dominated rate and such things must have been updated :), dominated rate scaled up for strictly positive!
		int dummIdx  		= 0;	// This will have to keep track of when real jumps are
		int dummIdx2		= 0;	// This will keep track of the frame
		double virtualRate 	= 0.0;
		double virtJumpTim	= 0.0; 
		boolean depNetwork	= true;
		conf dummyConf;

		// First configuration should be at time 0 and empty, only care about the frame!
		toReturn.nodes.add(new conf());
		toReturn.nodes.get(dummIdx).time = oldChain.get(dummIdx).time;
		toReturn.nodes.get(dummIdx).in = oldChain.get(dummIdx).in;
		toReturn.nodes.get(dummIdx).out = oldChain.get(dummIdx).out;

		outerloop1:
		while(true){

			// Work out rate for virtual jumps at this index, i.e. sum rates at all servers
			virtualRate = domRate * (1-arrivalsDomSum);
			for (int i=0; i<serverRates.size(); i++) {
				if (oldChain.get(dummIdx).configuration.get(i).size()>0) {
					virtualRate -= serverRates.get(i)[ oldChain.get(dummIdx).configuration.get(i).get(0).type ];
				}
			}
			dummyConf = oldChain.get(dummIdx); // Save current state to compare against actual transition

			// Now add the virtual jumps until next index time
			dummIdx++;
			dummIdx2++;
			ExponentialDistribution myExp = new ExponentialDistribution(rand2, 1/virtualRate);

			outerloop2:
			while(true){
				virtJumpTim = myExp.sample();
				if (virtJumpTim < oldChain.get(dummIdx).time - toReturn.nodes.get(dummIdx2-1).time) {
					toReturn.nodes.add(dummIdx2,new conf());
					toReturn.nodes.get(dummIdx2).time = toReturn.nodes.get(dummIdx2-1).time + virtJumpTim;
					dummIdx2++;
					// Add auxiliary variables 
					toReturn.aux.add( (rand2.nextDouble()>closeProb) ? new int[]  {-2,-1,-1} : new int[]  {-1,-1,-1} ); // For virtuals we use -1
				} else{
					break outerloop2;
				}
			}

			// Here actual arrival auxiliary variables
			if (oldChain.get(dummIdx).in != null) {
				for (int i=0; i<serverRates.size(); i++){
					if (oldChain.get(dummIdx).configuration.get(i).size()>0 &&
						oldChain.get(dummIdx).configuration.get(i).get(oldChain.get(dummIdx).configuration.get(i).size()-1).identifier
						== oldChain.get(dummIdx).in.identifier) {
						toReturn.nodes.add(dummIdx2,new conf());
						toReturn.nodes.get(dummIdx2).time = oldChain.get(dummIdx).time;
						toReturn.nodes.get(dummIdx2).in = oldChain.get(dummIdx).in;
						toReturn.aux.add( (rand2.nextDouble()>closeProb) ? new int[]  {-2,-1,-1} : new int[]  {oldChain.get(dummIdx).in.identifier,-1,i} );
					}
				}	
			} else { //Dep or transtion, focus on dummyConf servers (not empty, plus different to current server)
				for (int i=0; i<serverRates.size(); i++){
					if (dummyConf.configuration.get(i).size()>0 &&
						(oldChain.get(dummIdx).configuration.get(i).size() == 0 || 
							dummyConf.configuration.get(i).get(0).identifier != 
							oldChain.get(dummIdx).configuration.get(i).get(0).identifier)) {
						depNetwork = true;
						for (int j=0; j<serverRates.size(); j++){
							if (oldChain.get(dummIdx).configuration.get(j).size()>0 && 
								oldChain.get(dummIdx).configuration.get(j).get(oldChain.get(dummIdx).configuration.get(j).size()-1).identifier
								== dummyConf.configuration.get(i).get(0).identifier) {
								toReturn.nodes.add(dummIdx2,new conf());
								toReturn.nodes.get(dummIdx2).time = oldChain.get(dummIdx).time;
								toReturn.nodes.get(dummIdx2).served = oldChain.get(dummIdx).served;
								toReturn.nodes.get(dummIdx2).server = oldChain.get(dummIdx).server;
								toReturn.aux.add( (rand2.nextDouble()>closeProb) ? new int[]  {-2,-1,-1} : new int[]  {dummyConf.configuration.get(i).get(0).identifier,i,j} );
								depNetwork = false;
							}
						}
						if (depNetwork) { // Jump off the network
							toReturn.nodes.add(dummIdx2,new conf());
							toReturn.nodes.get(dummIdx2).time = oldChain.get(dummIdx).time;
							toReturn.nodes.get(dummIdx2).out = oldChain.get(dummIdx).out;
							toReturn.nodes.get(dummIdx2).server = oldChain.get(dummIdx).server;						
							toReturn.aux.add( (rand2.nextDouble()>closeProb) ? new int[]  {-2,-1,-1} : new int[]  {dummyConf.configuration.get(i).get(0).identifier,i,-1} );
						}
					}
				}
			}

			// Break main loop if we reach the end of chain
			if (dummIdx == oldChain.size() - 1)	break outerloop1;

		}
		
		return toReturn;

	}

	// Update a configuration not allowing arrivals or departures
	public void updateInnerConf(confProb dummyConfProb) {
		
		/*FCFS Policy*/
		// The departure we deal by rejection, just easier to code
		double prob = rand2.nextDouble()*(1-arrivalsDomSum);
		double sum 	= 0.0;

		for (int i=0; i<serverRates.size(); i++) {
			if (dummyConfProb.state.configuration.get(i).size() > 0){
				sum +=  serverRates.get(i)[dummyConfProb.state.configuration.get(i).get(0).type] / domRate;
				if (prob < sum) {
					//Transition from that server
					int destServer = probUtils.sampleIndex(transitions.get(dummyConfProb.state.configuration.get(i).get(0).type).get(i+1),rand2);
					if (destServer == serverRates.size()){ // NOOOOOOO Departure!!!! try again
						updateInnerConf(dummyConfProb);
						return;
					} else {
						dummyConfProb.upConnectProbs.add(transitions.get(dummyConfProb.state.configuration.get(i).get(0).type).get(i+1)[destServer] *
							serverRates.get(i)[dummyConfProb.state.configuration.get(i).get(0).type] / domRate);		
						dummyConfProb.state.configuration.get(destServer).add(dummyConfProb.state.configuration.get(i).get(0));
						dummyConfProb.state.configuration.get(i).remove(0);
						return;
					}
				}
			}
		}

		// If we reach here it means this is a virtual transition, yet we need the upConnectProbs
		dummyConfProb.upConnectProbs.add(1-arrivalsDomSum-sum);

	}



	/**
	Inner Classes for Jobs and System Configurations
	*/

	// A class to store a frame with jump times, entries/departures and auxiliaries
	class frame{

		public ArrayList<conf> nodes 	= new ArrayList<conf>();
		public ArrayList<int[]> aux	= new ArrayList<int[]>(); // So... identifier, dep queue, entry queue

	}

	// A class to keep track of probabilities linked to states in Filtering procedure
	class confProb {

		public conf state;
		public double stateProb;
		public ArrayList<Integer> upConnect; 		// Connectivity to previous layer in Filter
		public ArrayList<Double> upConnectProbs; 	// Connectivity to previous layer in Filter, Probabilities

		public confProb(){
			upConnect 		= new ArrayList<Integer>();
			upConnectProbs 	= new ArrayList<Double>();
		}

	}

	// A class of Observations in data
	class obs {
		public double 				time 		= 0;
		public int 					in 			= -1;
		public int 					inType 		= -1;		
		public int 					out 		= -1;
		public int 					served 		= -1;
		public int 					server 		= -1;		
		public ArrayList<Integer> 	jobsInIDs 	= new ArrayList<Integer>();
	}


	// A class of Jobs, including a unique identifier and its type
	class job {

		public int identifier;
		public int type;
		public job (int id, int typ){
			identifier = id;
			type = typ;
		}

	}

	// A class of objects <time, system configuration>
	class conf implements Cloneable {

		public double time;
		public job in;  // To keep record of when jobs come to system
		public job out;	// To keep record of when jobs leave the system
		public int served = -1;  // To keep record of served jobs, at whatever server
		public int server = -1;	// To keep record of servers at which jobs are served
		public ArrayList<ArrayList<job>> configuration = new ArrayList<ArrayList<job>>(serverRates.size());
		public conf clone() throws CloneNotSupportedException {
			return (conf) super.clone();
		}
		public conf(){
			for (int i=0; i<serverRates.size(); i++) configuration.add(new ArrayList<queueNetwork.job>());
		}

	}

	/**
	Constructor
	*/

	// Define a Constructor Providing Queue Form and Parameters
	public queueNetwork(int jobTypes, int servers, Double[] arr , Double[][] rat, Double[][][] trans){

		// Add arrival rates
		for (int i=0; i<jobTypes; i++ ) {
			arrivals.add(arr[i]);
		}
		// Add server Rates
		for (int i=0; i<servers; i++ ) {
			serverRates.add(rat[i]);
		}
		// Add transition Matrices
		for (int i=0; i<jobTypes; i++) {
			ArrayList<Double[]> dummy = new ArrayList<Double[]>();
			for (int j=0; j<(servers+1); j++) {
				dummy.add(trans[i][j]);
			}
			transitions.add(dummy);
		}
		// Work out Dominating Rate for Uniformization
		for (Double i : arr){
			domRate += i;
		}
		for (Double[] i : rat){
			domRate += Collections.max(Arrays.asList(i));
		}
		// Scale it up
		domRate *= scaleDomRate;
		// Add dominated arrival rates
		for (int i=0; i<jobTypes; i++ ) {
			arrivalsDom.add(arr[i]/domRate);
		}
		System.out.println("\nNetwork Updated, the strictly greater dominating rate is: " + domRate + "\n");
		//Sum them over to avoid cumbersome operations later
		for (int j=0; j<arrivalsDom.size(); j++) {
			arrivalsDomSum += arrivalsDom.get(j);
		}

	}
}
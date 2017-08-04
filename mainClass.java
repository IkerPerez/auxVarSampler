
import java.util.ArrayList;
import java.util.List;
import java.io.FileWriter;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.stream.Collectors;

public class mainClass {
		   
	public static void main (String args[]) {
							   
		queueNetwork bottleneckQueue = new queueNetwork(
		3, 									// Different Job Types in system
		3, 									// Number of Servers in System
		new Double[] {0.08,0.06,0.04}, 		// Instantiate arbitrary Arrival Rates
		new Double[][] {{0.3,	0.25, 0.2},
						{0.7,	0.5, 0.3}, 
						{1.5,	1.2, 0.8}},	// Instantiate arbitrary Server Rates for Job types (respect prior inequalities)
		new Double[][][] { {{0.5	,0.5	,0.0	,0.0},
			 				{0.0	,0.0	,1.0	,0.0},
			 				{0.0	,0.0	,1.0	,0.0},
			 				{0.0	,0.0	,0.0	,1.0}}, // Instantiate Job Type 1 Transitions 
						   {{0.5	,0.5	,0.0	,0.0},
			 				{0.0	,0.0	,1.0	,0.0},
			 				{0.0	,0.0	,1.0	,0.0},
			 				{0.0	,0.0	,0.0	,1.0}}, // Instantiate Job Type 2 Transitions
						   {{0.5	,0.5	,0.0	,0.0},
			 				{0.0	,0.0	,1.0	,0.0},
			 				{0.0	,0.0	,1.0	,0.0},
			 				{0.0	,0.0	,0.0	,1.0}}	// Instantiate Job Type 3 Transitions
						 }
		);

	 	/**
	 	Write and output a small set of simulations from the above queue; 5 arbitrary units of time.
		*/

 		// for (int i = 0; i<500 ; i++) {
 		// 	bottleneckQueue.writeSimulationUnbounded(5,"\\Files\\" + "simulation_" + String.valueOf(i+1) + ".csv");
 		// }

 		/**
 		Produce output file assuming full observations, por testing purposes if needed. Seeds are fixed, to match output above comment the previous out.
 		*/

    	ArrayList<ArrayList<queueNetwork.conf>> iepe = new ArrayList<ArrayList<queueNetwork.conf>>();
    	for (int i=0; i<500; i++) iepe.add(bottleneckQueue.doSimulation(5));
 		try{
 			FileWriter writer = new FileWriter(System.getProperty("user.dir") + "\\OutputReal.csv");
 			writer.append("aj1, aj2, aj3, s1j1, s1j2, s1j3, s2j1, s2j2, s2j3, s3j1, s3j2, s3j3, \n");
 			for (int jj = 0; jj<1000 ; jj++) {
 				System.out.println( "Progress:" + jj/1000.0 * 100.00 + "%");

				// Update
 				bottleneckQueue.updateRates(iepe);
 				
				// Write to file
				for (int i=0; i<bottleneckQueue.arrivals.size(); i++) writer.append(bottleneckQueue.arrivals.get(i) + ", ");
					for (int serv=0; serv<bottleneckQueue.serverRates.size(); serv++) {
						for (int i=0; i<bottleneckQueue.arrivals.size(); i++) {
							writer.append(bottleneckQueue.serverRates.get(serv)[i] + ", "); //Screwed due to first sequences forced by me :)
						}
					}
				writer.append("\n");
			}
			writer.flush();
			writer.close();
		}catch (IOException e){
			e.printStackTrace();
		}

		/**
		Sampler
 		*/

		// Read and manipulate the simulations for the purposes of drawing inference on input parameters
		ArrayList<ArrayList<queueNetwork.obs>> simulations = new ArrayList<ArrayList<queueNetwork.obs>>();
		for (int i = 0; i<500 ; i++) {
			simulations.add(bottleneckQueue.readSimulation("simulation_" + String.valueOf(i+1) + ".csv"));
		}

		// Initial arbitrary network path sequences given observation
		ArrayList<ArrayList<queueNetwork.conf>> stSeqs = new ArrayList<ArrayList<queueNetwork.conf>>();
		stSeqs = simulations.parallelStream()		
			.map(i  -> bottleneckQueue.initFFBS(i) ) 	
			.collect(Collectors.toCollection(ArrayList::new));

		//Adapt the initialized paths to your starting values, to reduce the initialization bias
		for(int j=0; j<1500; j++){
			System.out.println(j);
			stSeqs = stSeqs.parallelStream()		
				.map(i  -> bottleneckQueue.newChain(i) ) 	
				.collect(Collectors.toCollection(ArrayList::new));	
		}

		//Print current rate estimates, will rarely match initial values since may be infeasible given network topology
		System.out.println("Starting Sequences adapted to starting values");
		bottleneckQueue.updateRates(stSeqs);
		System.out.print("\n");
		for (int i=0; i<bottleneckQueue.arrivals.size(); i++) System.out.print(bottleneckQueue.arrivals.get(i) + "\t");
			System.out.print("\n\n");
		for (int serv=0; serv<bottleneckQueue.serverRates.size(); serv++) {
			for (int i=0; i<bottleneckQueue.arrivals.size(); i++) {
				System.out.print(bottleneckQueue.serverRates.get(serv)[i] + "\t"); 
			}
			System.out.print("\n");
		}

		// MCMC sampler
		try{
			FileWriter writer = new FileWriter(System.getProperty("user.dir") + "\\Output.csv");
			writer.append("aj1, aj2, aj3, s1j1, s1j2, s1j3, s2j1, s2j2, s2j3, s3j1, s3j2, s3j3, \n");
			int itersMCMC = 0;
			for (int jj = 0; jj<100000 ; jj++) {
				// Progress info:
				System.out.println( "Progress:" + jj/100000.0 * 100.00 + "%");
				// New network paths
    			stSeqs = stSeqs.parallelStream()		
    				.map(i  -> bottleneckQueue.newChain(i) ) 	
    				.collect(Collectors.toCollection(ArrayList::new));	
				// Update Rates
				bottleneckQueue.updateRates(stSeqs);
				// Write to file
				for (int i=0; i<bottleneckQueue.arrivals.size(); i++) writer.append(bottleneckQueue.arrivals.get(i) + ", ");
					for (int serv=0; serv<bottleneckQueue.serverRates.size(); serv++) {
						for (int i=0; i<bottleneckQueue.arrivals.size(); i++) {
						writer.append(bottleneckQueue.serverRates.get(serv)[i] + ", "); 
					}
				}
				writer.append("\n");
			}
			writer.flush();
			writer.close();
		}catch (IOException e){
			e.printStackTrace();
		}
 		
	}
}
 
 


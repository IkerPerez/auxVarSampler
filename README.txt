

### INTRODUCTION ###

This repository includes data and code allowing to reproduce results within the following manuscript:

Perez, I., Hodge, D. and Kypraios, T. – Auxiliary variables for Bayesian inference in multi-class queueing networks. https://arxiv.org/abs/1703.03475.


### SAMPLE DATA ###

Within the "Files" folder, we include 500 sample simulations from a bottleneck network as introduced in Example 2 in the paper. Real data can replace these files, shold the structure remain the same.


### JAVA CODE ###

queueNetwork.java includes methods and instance variables that allow to simulate data or calibrate the MCMC sampler. Of relevance:

- Instance variables:

1) scaleDomRate: Allows to determine the dominating rate multiplier in the sampler.
2) closeProb2: Determines the likelihood for "campling" nodes within the sampler.
3) fixedTrans: Used when we want to draw inference on network topology structure (see updateRates below).
4) missingData: determines the approximate percentage of data within the simulations that we will consider "missing".

- Methods:

1) newChain: taking a prior network path as an origin, it will conduct a full iteration in the sampler and produce a new path compatible with the observations.
2) updateRates: given a collection of network path sequences, it will produce sampler service rates across the stations. The prior specification must be included within this method (see the paper). If we wish to draw inference on network topology, lines must be commented out and manipulated accordingly.
3) initFFBS: method to instantiate arbitrary network paths given observations; in order to initialize MCMC sampler.
4) readSimulation: allows to read the simulations; in addition, it may read any source of real data that resembles the structure of the simulations included in the "Files" folder.
5) writeSimulationUnbounded: allows to simulate network paths and observations from a user-specified network.

mainClass.java is the executable file. The network structure and initialization parameters must be fed to this file in a way that resembles the current version of the code.


For any queries, please contact Dr. Iker Perez at iker.perez@nottingham.ac.uk.
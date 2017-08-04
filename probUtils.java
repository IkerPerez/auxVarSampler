
import org.apache.commons.math3.random.SynchronizedRandomGenerator;

public class probUtils {

	public static int sampleIndex(Double[] vector, SynchronizedRandomGenerator rand){

		double p 	= rand.nextDouble();
		double sum	= 0;
		for (int i=0; i<vector.length; i++) {
			sum += vector[i];
			if (p<sum) return i;
		}
		return vector.length+1;
	}

}

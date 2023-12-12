package section6_1;
import java.util.ArrayList;
import java.util.Random;

public class main {
	public static void main(String[] args) {
		int repeat = 1000;
		int correct1 = 0;
		int correct2 =0;
		int correct3 =0;
		double alpha = 0.05;
		double sampleSize1 = 0d;
		double sampleSize2 = 0d;
		double sampleSize3 = 0d;
		long seed= 100;
		int k = 20;
		int n0 = 10;
		double delta = 2d;
		Random R = new Random(seed);
		
		
		
		
		
		for(int count = 0; count < repeat; count++) {
			ArrayList<Double> mu = new ArrayList<Double>();
			ArrayList<Double> sigma2 = new ArrayList<Double>();
			for(int i = 0; i < k; i++) {
				if(i!=0) {
					//mu.add(0d);
					mu.add(R.nextDouble()*0.1*k);
					//mu.add(1.5-0.5*(i+1));
				}else {
					mu.add(0.1*k+0.5);
					//mu.add(1d);
				}
				//mu.add(0.4*i);
				//sigma2.add(2d);
				sigma2.add(R.nextDouble()*10+5);
				//sigma2.add(10d);
			}
			
			//mu.add(0.01d);mu.add(0d);
			//sigma2.add(1d);sigma2.add(1d);
			
			long tempSeed = R.nextLong();
			//FHN1 y1 = new FHN1(10,alpha,mu,sigma2,tempSeed);//R.nextLong());
			procedure2 y1 = new procedure2(alpha,mu,sigma2,tempSeed);//R.nextLong());
			//procedure1 y1 = new procedure1(alpha,mu,sigma2,tempSeed);
			//KN y1 = new KN(n0,delta,alpha,mu,sigma2,tempSeed);
			y1.run();
			//y2.run();
			//y3.run();
			if(y1.getBestID()==0) {
				correct1 +=1;
			}
			//if(y2.getBestID()==0) {
				//correct2 +=1;
			//}
			//if(y3.getBestID()==0) {
				//correct3 +=1;
			//}
			sampleSize1 += y1.getTotalSampleSize();
			//sampleSize2 += y2.getTotalSampleSize();
			//sampleSize3 += y3.getTotalSampleSize();
			//System.out.println(count+" "+y1.getBestID()+" "+y1.getTotalSampleSize()+" "+y2.getBestID()+" "+y2.getTotalSampleSize());
			System.out.println(y1.getBestID()+" "+y1.getTotalSampleSize());
			//System.out.println(count+" "+y3.getBestID()+" "+y3.getTotalSampleSize());	
		}
		
		System.out.println(correct1*1d/repeat+" "+sampleSize1/repeat);
		//System.out.println(correct2*1d/repeat+" "+sampleSize2/repeat);
		//System.out.println(correct3*1d/repeat+" "+sampleSize3/repeat);
	}
}

package section6_2;

import java.util.ArrayList;
import java.util.Random;

import org.apache.commons.math3.distribution.TDistribution;

public class P2 {
	double alpha,delta;
	int bestID=-1,sampleSize=0;
	int seed = 1374;
	double bestAVG = 0d;
	int bestSampleSize = 0;
	int sp,ep;
	double muCall, muServe;
	int k;
	
	int goodID = -1;
	
	double[][] altVec;
	
	Random R = new Random();
	
	public P2() {
		alpha = 0.05;
		delta = 0.1;
		seed = 123;
		sp = 1000;
		ep = 2000;
		muCall = 55d/60;
		muServe = 10d/60;
	}
	public P2(int sp, int ep,double muCall, double muServe,
			double alpha, double delta,double[][] altVec, int seed) {
		this.sp = sp;
		this.ep = ep;
		this.muCall = muCall;
		this.muServe = muServe;
		this.alpha = alpha;
		this.delta = delta;
		this.altVec = altVec;
		this.seed = seed;
		
	}
	
	//Return best alternative
	public int getBestID() {
		return bestID;
	}
	
	
	//Return total sample size (include the sample size of the reference alternative)
	public int getSampleSize() {
		return sampleSize;
	}
	
	
	//Return the sample average of the best alternative
	public double getBestAVG(){
		return bestAVG;
	}
	
	//Return the sample size of the best alternative
	public int getBestSampleSize() {
		return bestSampleSize;
	}
	
	
	
	//While conducting the selection, the reference alternative is the 1st alternative in set I.
	public void run() {
		double[] tempAlt = new double[9];
		R.setSeed((long)seed);
		genObv genObv = new genObv(sp,ep,muCall,muServe,R.nextInt());
		ArrayList<Integer> I = new ArrayList<Integer>();
		k = altVec.length;
		for(int i = 0 ; i < k; i++) {
			I.add(i);
		}
		 double[][] Xij =  new double[k][k];
		 double[][] Xij2 = new double[k][k];
		 int r  = 1;
		 while(I.size()>1) {
			 int t = (int)Math.pow(2, r);
			 
			 int n =  t/2;
			 
			 if(r == 1) {
				 n = 2;
			 }
			 
			 
			 double[] tempX = new double[I.size()];
			 
			 for(int ell = 0; ell < n; ell++) {
				 for(int i = 0; i < I.size(); i++) {
					 for(int _e=0; _e < 9; _e++) {
							tempAlt[_e] =altVec[I.get(i)][_e];
						}
						
					genObv.setAlt(tempAlt);
					genObv.run();
					//if(I.get(i)==54) {
						
						//System.out.println(genObv.getPerOfGC());
					//}
					/**if(I.get(i)==35) {
						
						System.out.println(genObv.getPerOfGC());
					}**/
					tempX[i]=genObv.getPerOfGC();
				 }
				 
				 for(int i = 0; i < I.size(); i++) {
					 Xij[I.get(i)][I.get(i)] = Xij[I.get(i)][I.get(i)] + tempX[i];
					 Xij2[I.get(i)][I.get(i)] = Xij2[I.get(i)][I.get(i)] + tempX[i] * tempX[i];
					 
					 for(int j = i + 1; j < I.size(); j++) {
						 Xij[I.get(i)][I.get(j)] = Xij[I.get(i)][I.get(j)]+tempX[i]-tempX[j];
						 Xij[I.get(j)][I.get(i)] = Xij[I.get(j)][I.get(i)]+tempX[j]-tempX[i];
						 Xij2[I.get(i)][I.get(j)] = Xij2[I.get(i)][I.get(j)] +(tempX[i]-tempX[j])*(tempX[i]-tempX[j]);
						 Xij2[I.get(j)][I.get(i)] = Xij2[I.get(i)][I.get(j)];
					 }
				 }
			 }
		
			 for(int i  = 0;  i < I.size(); i++) {
				 boolean check  = true;
				 for (int j = 0; j < I.size(); j++) {
					 if(j!=i) {
						 double Sij2 = (Xij2[I.get(i)][I.get(j)] -  Xij[I.get(i)][I.get(j)]* Xij[I.get(i)][I.get(j)]/t)/(t-1);
						 double Zij = Xij[I.get(i)][I.get(j)]/Math.sqrt(Sij2);
						 if(Zij+delta*t/Math.sqrt(Sij2)<hht(t)) {
							 check = false;
							 break;
						 }
					 }
					 
				 }
				 if(check) {
					 double tempLargest = 0d;
					 int index = 0;
					 
					 for(int _s = 1; _s<I.size();_s++) {
						// System.out.println(I.get(0)+" "+I.get(_s)+" "+Xij[I.get(_s)][I.get(0)]);
						 if(Xij[I.get(_s)][I.get(0)]>tempLargest) {
							 index = _s;
							 tempLargest = Xij[I.get(_s)][I.get(0)];
						 }
					 }
					 
					 
					 sampleSize = sampleSize + t * (I.size()-1);
					 bestID = I.get(index);
					 I.clear();
					 I.add(bestID);
					 break;
				 }
			 } 
				
			
			 
			 for(int i  = 0;  i < I.size(); i++) {
				 for (int j = 0; j < I.size(); j++) {
					 if(i!=j) {
						 double Sij2 = (Xij2[I.get(i)][I.get(j)] -  Xij[I.get(i)][I.get(j)]* Xij[I.get(i)][I.get(j)]/t)/(t-1);
						 double Zij = Xij[I.get(i)][I.get(j)]/Math.sqrt(Sij2);
						/** if(I.get(j)==54&&I.get(i)==35) {
							 System.out.println(Zij+" "+Xij2[I.get(i)][I.get(j)]+" "+Xij[I.get(i)][I.get(j)]* Xij[I.get(i)][I.get(j)]/t);
						 }**/
						 //System.out.println(t+" "+Zij+" "+ht(t)+" "+Math.sqrt(Sij2)+" "+Xij[I.get(i)][I.get(j)]);
						 if(Zij > ht(t)&&Sij2!=0) {
							 
							 sampleSize = sampleSize + t;
							 I.remove(j);
							 if(j < i) {
								 i--;
							 }
							 j--;
						 }
					 }
				 }
			 }
			 if(I.size()==1) {
				 sampleSize = sampleSize + t;
				// System.out.println(t);
			 }
			/**System.out.println(" ");
			 for(int i  = 0;  i < I.size(); i++) {
				 System.out.print(I.get(i)+" ");
			 }
			 System.out.println(" ");**/
			 r++;
		 }
		 bestID = I.get(0);
		
		
		
		
		
	
	}
	public double hht(int t) {
		double beta = alpha*(2-Math.PI*Math.PI/6)/(Math.PI*Math.PI/6-1);
		double p =  1-beta/((k-1)*Math.pow(Math.log(2*t)/Math.log(2),2));
		 TDistribution td = new TDistribution(t-1);
		 
		 double  hht = Math.sqrt(t)*td.inverseCumulativeProbability(p);
		 return hht;
	}
	
	public double ht(int t) {
		
		 
		 double p =  1-alpha/((k-1)*Math.pow(Math.log(2*t)/Math.log(2),2));
		 TDistribution td = new TDistribution(t-1);
		 
		 double  ht = Math.sqrt(t)*td.inverseCumulativeProbability(p);
		 
		 return ht;
	 }
}

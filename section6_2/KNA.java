package section6_2;

import java.util.ArrayList;
import java.util.Random;

public class KNA {
	double alpha,delta;
	int bestID=-1,sampleSize=0;
	int seed = 1374;
	double bestAVG = 0d;
	int bestSampleSize = 0;
	int n0Prime;
	int sp,ep;
	double muCall, muServe;
	
	double[][] altVec;
	
	Random R = new Random();
	
	public KNA() {
		alpha = 0.05;
		delta = 0.1;
		seed = 123;
		n0Prime = 10;
		sp = 1000;
		ep = 2000;
		muCall = 55d/60;
		muServe = 10d/60;
	}
	public KNA(int sp, int ep,double muCall, double muServe,
			double alpha, double delta,double[][] altVec, int seed, 
			int n0Prime) {
		this.sp = sp;
		this.ep = ep;
		this.muCall = muCall;
		this.muServe = muServe;
		this.alpha = alpha;
		this.delta = delta;
		this.altVec = altVec;
		this.seed = seed;
		this.n0Prime  = n0Prime;
		
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
		int k = altVec.length;
		for(int i = 0 ; i < k; i++) {
			I.add(i);
		}
		double[] X = new double[I.size()];
		
		double[][] tempX = new double[I.size()][n0Prime];
		double[][] S2 = new double[I.size()][I.size()];
		

		double hr2 = (n0Prime - 1)*(Math.pow(2*alpha/(k-1), -2.0/(n0Prime-1))-1)/(2*delta);
		for(int t = 0;  t < n0Prime;t++) {
			for(int i=0; i < I.size(); i++) {
				
				for(int ell=0; ell < 9; ell++) {
					tempAlt[ell] =altVec[i][ell];
				}
					
				genObv.setAlt(tempAlt);
				genObv.run();
					
				tempX[i][t]=genObv.getPerOfGC();
				X[i] = X[i] + tempX[i][t];
					
				
			}
		}
		
		for(int i = 0 ; i <I.size(); i++) {
			for(int j = i + 1; j < I.size(); j++) {
				double tempS2 = 0;
				double tempAVGDiff = X[i]/n0Prime - X[j]/n0Prime;
				for(int count = 0 ; count < n0Prime; count++) {
					tempS2 =  tempS2 + (tempX[i][count]-tempX[j][count]-tempAVGDiff)* (tempX[i][count]-tempX[j][count]-tempAVGDiff);
				}
				tempS2 = tempS2 / (n0Prime - 1.0);
				S2[i][j]=tempS2;
				S2[j][i]=tempS2;
			}
		}
		int t = n0Prime;
		while(I.size()>1) {
			for(int i = 0 ; i < I.size(); i++) {
				boolean checkEli = false;
				//if(I.get(i)==20) {
				//	System.out.println("Sample Average: "+I.get(i)+" "+X[I.get(i)]/t);
				//}
				for(int j = i + 1 ; j < I.size(); j++) {
					
					double boundary=hr2*S2[I.get(i)][I.get(j)] - delta * t/2;
					double Zij = X[I.get(i)]-X[I.get(j)];
					if((boundary > 0&&Zij>boundary)||(boundary<0&&Zij>0)) {
						
						
						sampleSize =sampleSize + t;
						
						I.remove(j);						
						j--;
						
					}
					
					if((boundary > 0&&Zij<-boundary)||(boundary<0&&Zij<0)) {
						
						
						checkEli = true;
						
					}
				}
				
				if(checkEli && I.size()>1) {
					
					
					sampleSize = sampleSize + t;
					
					I.remove(i);
				}
			}
			for(int i = 0 ; i < I.size(); i++) {
				
					for(int ell=0; ell < 9; ell++) {
						tempAlt[ell] =altVec[I.get(i)][ell];
					}
					
					genObv.setAlt(tempAlt);
					genObv.run();
					
					X[I.get(i)] = X[I.get(i)] + genObv.getPerOfGC();
				
			}
			//for(int _s = 0; _s <I.size();_s++) {
			//	System.out.print(I.get(_s)+" ");
			//}
			t++;
		}
		sampleSize = sampleSize + t;

		
		bestAVG = X[I.get(0)]/t;
		bestSampleSize = t;
		bestID = I.get(0);
	}
	
}

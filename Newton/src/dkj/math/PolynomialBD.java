package dkj.math;

import java.util.ArrayList;
import java.math.BigDecimal;
import dkj.math.NoIVException;

/**************************************************************************
 * 
 * PolynomialBD  -  This class exemplifies a mathematical polynomial of
 * 	the form a0 + a1*x + a2*x^2 + ... .  It uses an ArrayList to store
 * 	the coefficients in increasing order of power until the highest
 *  non-zero term.  It provides methods for adding, multiplying, taking
 *  derivatives and finding the value of the polynomial at a given point.
 *  It also computes the degree of the polynomial, allows the modification
 *  of individual terms and provides a means to obtain a human-readable
 *  string. The variable can be set for display purposes only.  Changing
 *  the variable has NO EFFECT on multiplication or addition.
 *
 *************************************************************************/

public class PolynomialBD {
	private ArrayList<BigDecimal> alCoeffs;
	private String sVar = "x";
	private static int scale = 32;
	private static boolean DEBUG = false;
	
	public PolynomialBD() {
		alCoeffs = new ArrayList<BigDecimal>();
	}
	
	public PolynomialBD(BigDecimal[] dCoeffs) {
		alCoeffs = new ArrayList<BigDecimal>();
		for (int i=0; i<dCoeffs.length; i++)
			alCoeffs.add(dCoeffs[i]);
	}
	
	public PolynomialBD(int iDegree) {
		alCoeffs = new ArrayList<BigDecimal>();
		for (int i=0; i<=iDegree; i++)
			alCoeffs.add(BigDecimal.valueOf(0.0));
	}
	
	public static void setScale(int iScale) {
		scale = iScale;
	}

	public void setTerm(int degree, BigDecimal coeff) {
		if (alCoeffs.size()<(degree+1)) {
			while (alCoeffs.size()<degree)
				alCoeffs.add(BigDecimal.valueOf(0.0));
		alCoeffs.add(coeff);
		} else {
			alCoeffs.set(degree,coeff);
		}
	}
	
	public int degree() {
		return alCoeffs.size()-1;
	}
	
	public void add(PolynomialBD pSecond){
		if (this.degree()<pSecond.degree()) {
			for (int i=0; i<pSecond.degree()-this.degree(); i++)
				alCoeffs.add(BigDecimal.valueOf(0.0));
		}
		int j=0;
		for (BigDecimal i : pSecond.alCoeffs) {
			this.alCoeffs.set(j, this.alCoeffs.get(j).add(i));
			j++;
		}
	}
	
	public static PolynomialBD add(PolynomialBD pFirst, PolynomialBD pSecond) {
		pFirst.add(pSecond);
		return pFirst;
	}
	
	public void multiply(PolynomialBD pSecond) {
		PolynomialBD pTemp =new PolynomialBD(this.degree()+pSecond.degree());
		for (int i=0; i<=this.degree(); i++)
			for (int j=0; j<=pSecond.degree(); j++) {
				pTemp.alCoeffs.set(i+j, 
					pTemp.alCoeffs.get(i+j).add(this.alCoeffs.get(i).multiply(pSecond.alCoeffs.get(j))));
			}
		this.set(pTemp);
	}
	
	public static PolynomialBD multiply(PolynomialBD pFirst, PolynomialBD pSecond) {
		pFirst.multiply(pSecond);
		return pFirst;
	}
	
	public void set(PolynomialBD pSecond) {
		this.clear();
		for (BigDecimal i: pSecond.alCoeffs)
			this.alCoeffs.add(i);
	}
	
	public void clear() {
		while (!this.alCoeffs.isEmpty())
			this.alCoeffs.remove(0);
	}

	public void differentiate() {
		for (int i=0; i<this.degree(); i++)
			this.alCoeffs.set(i, this.alCoeffs.get(i+1).multiply(BigDecimal.valueOf(i+1)));
		this.alCoeffs.remove(this.degree());
	}
	
	public static PolynomialBD differentiate(PolynomialBD pFirst) {
		PolynomialBD pTemp = new PolynomialBD();
		pTemp.set(pFirst);
		pTemp.differentiate();
		return pTemp;
	}
	
	public BigDecimal valueAt(BigDecimal x) {
		BigDecimal y=BigDecimal.valueOf(0);
		for (int i=0; i<=this.degree(); i++)
			y= y.add(this.alCoeffs.get(i).multiply(x.pow(i)));
		return y;
	}
	
	public PolynomialBD tangentLineAt(BigDecimal x) {
		PolynomialBD pTemp=new PolynomialBD();
		BigDecimal y = this.valueAt(x);
		PolynomialBD pDeriv = PolynomialBD.differentiate(this);
		BigDecimal m = pDeriv.valueAt(x);
		pTemp.setTerm(0, y.subtract(m.multiply(x)));
		pTemp.setTerm(1, m);
		return pTemp;
	}
	
	public BigDecimal tangentRootNear(BigDecimal x) {
		PolynomialBD pDeriv = PolynomialBD.differentiate(this);
		BigDecimal y = this.valueAt(x);
		BigDecimal yprime = pDeriv.valueAt(x);
		return (yprime.multiply(x).subtract(y)).divide(yprime,scale,BigDecimal.ROUND_HALF_EVEN);		
	}
	
	public BigDecimal bisectionRoot (BigDecimal bdLower, BigDecimal bdUpper, BigDecimal bdEpsilon) throws NoIVException {
		BigDecimal bdYLower, bdYUpper, bdMidpoint, bdYMidpoint;
		int i=0;
		bdYLower=this.valueAt(bdLower);
		bdYUpper=this.valueAt(bdUpper);
		if (bdYUpper.multiply(bdYLower).compareTo(BigDecimal.ZERO)>-1) {
			throw new NoIVException();
		} else {
			do {
				i++;
				bdMidpoint = bdLower.add(bdUpper).divide(new BigDecimal(2.0));
				bdYMidpoint = this.valueAt(bdMidpoint);
				if (bdYMidpoint.abs().compareTo(bdEpsilon)<0) {
					bdLower=bdMidpoint;
					bdUpper=bdMidpoint;
				} else if (bdYUpper.multiply(bdYMidpoint).compareTo(BigDecimal.ZERO)>-1) {
					assert (bdYLower.multiply(bdYMidpoint).compareTo(BigDecimal.ZERO)<0);
					bdUpper=bdMidpoint;
				} else {
					assert (bdYUpper.multiply(bdYMidpoint).compareTo(BigDecimal.ZERO)<0);
					bdLower=bdMidpoint;
				}
				if (DEBUG)
					System.out.println("Iteration: "+i+" Lower: "+bdLower+" Upper: "+bdUpper);
			} while(bdUpper.subtract(bdLower).compareTo(bdEpsilon)>-1);
		}
		return bdLower;
	}
	
	public void setVar(String newVar) {
		sVar = newVar;
	}
	
	public String toString() {
		String poly="";
		int exp=0;
		for (BigDecimal i : alCoeffs) {
			if (!i.equals(BigDecimal.valueOf(0.0))) {
				poly = i.toString() + getVar(sVar,exp) 
					+ getSign(poly,exp) + poly;
			}
			exp++;
		}
		return poly;
	}
	
	String getVar(String var, int exp) {
		if (exp==0)
			return "";
		else if (exp==1)
			return var;
		else
			return var+"^"+exp;
	}
	
	String getSign(String poly, int exp) {
		if (poly.isEmpty() || poly.startsWith("-"))
			return "";
		else
			return "+";
	}
}

package newton;

import java.math.BigDecimal;
import dkj.math.PolynomialBD;
import dkj.math.NoIVException;

public class Newton {
	
	static BigDecimal epsilon = new BigDecimal(1e-293);

	public static void main(String[] args) {
		PolynomialBD.setScale(300);
		PolynomialBD pMine =  new PolynomialBD(new BigDecimal[] 
				{BigDecimal.valueOf(-7),BigDecimal.valueOf(0),BigDecimal.valueOf(1)});
		BigDecimal dGuess=BigDecimal.valueOf(2);
		int i=0;
		System.out.println("Iter: "+i+"  Guess: "+dGuess+"  Value: "+pMine.valueAt(dGuess));
		while(pMine.valueAt(dGuess).abs().compareTo(epsilon)>-1) {
			i++;
			dGuess=pMine.tangentRootNear(dGuess);
			System.out.println("Iter: "+i+"  Guess: "+dGuess+"  Value: "+pMine.valueAt(dGuess));
		}
		try {
			System.out.println("Zero: " + pMine.bisectionRoot(
				BigDecimal.valueOf(2), BigDecimal.valueOf(3), epsilon));
		} catch (NoIVException e){
			System.out.println("We threw an exception!");
		}
	}

}

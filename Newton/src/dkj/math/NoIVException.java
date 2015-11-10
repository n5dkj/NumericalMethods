package dkj.math;

public class NoIVException extends Exception {
	static final long serialVersionUID = 31415925;
	
	public NoIVException() {
		super("There is no sign change in the range and the intermediate value theorem does not apply!");
	}
}

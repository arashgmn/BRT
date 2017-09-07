package objects;

import java.util.ArrayList;

public interface Fillable
{
	public static ArrayList<Fillable> all = new ArrayList<Fillable>();
	public ArrayList<Person> passengers = new ArrayList<Person>();
	public int maxCapacity();
	public int occupiedCapacity();
	public int remainingCapacity();
	public int type();
}

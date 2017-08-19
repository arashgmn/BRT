package objects;

import java.util.ArrayList;

public interface Fillable
{
	public ArrayList<Person> passengers = new ArrayList<Person>();
	public int maxCapacity();
	public int occupiedCapacity();
	public int remainingCapacity();
	void addPerson();
}

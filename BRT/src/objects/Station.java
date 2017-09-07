package objects;

public class Station implements Fillable //type 2
{
	private int capacity=150;
	private double position=0;
	
	public Station(int capacity, double position)
	{
		this.position=position;
		this.capacity=capacity;
		Fillable.all.add(this);
	}
	
	public double position()
	{
		return position;
	}

	@Override
	public int maxCapacity()
	{
		return capacity;
	}

	@Override
	public int occupiedCapacity()
	{
		return passengers.size();
	}

	@Override
	public int remainingCapacity()
	{
		return capacity-passengers.size();
	}

	@Override
	public int type()
	{
		return 2;
	}
}

package objects;

public class Bus implements Fillable //type 1
{
	private int capacity=40;
	private double position=0;
	private double speed=10;
	public Station station=null;

	public Bus(int capacity, double position)
	{
		this.capacity=capacity;
		this.position=position;
		Fillable.all.add(this);
	}
	
	public double getSpeed()
	{
		return speed;
	}
	
	public void setSpeed(double s)
	{
		this.speed=s;
	}
	
	public double getPosition()
	{
		return position;
	}
	
	public void setPosition(double p)
	{
		this.position=p;
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
		return 1;
	}
}

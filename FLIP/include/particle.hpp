
#define NUM_PARTICLES 4
#define SQRT_NUM_PARTICLES 2

struct Point
{
    Point(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    };

    Point()
    {
        this->x = 0;
        this->y = 0;
        this->z = 0;
    };

    double x;
    double y;
    double z;

    double l1_distance(Point p1)
    {
        return abs(x - p1.x) + abs(y - p1.y) + abs(z - p1.z);
    }

    double bilinear(Point p1)
    {
        return abs(x - p1.x)*abs(y - p1.y);
    }
};

struct Particle
{
    public:
        Particle(int id, Point position, Point velocity, double mass)
        {
            // TODO make sure this is memory safe in terms of ownership
            this->velocity = velocity;
            this->position = position;
            this->id = id;
            this->mass = mass;
        };

        Particle(int id, double mass)
        {
            // TODO make sure this is memory safe
            this->velocity = Point();
            this->position = Point();
            this->id = id;
            this->mass = mass;
        };

        Particle()
        {
            // TODO make sure this is memory safe
            this->velocity = Point();
            this->position = Point();
            this->id = 0;
            this->mass = 0;
        };

        void advect(double timestep)
        {
            this->position.x += timestep*velocity.x;
            this->position.y += timestep*velocity.y;
            this->position.z += timestep*velocity.z;
        }

        void apply_force(Point force, double timestep)
        {
            this->velocity.x += timestep*force.x;
            this->velocity.y += timestep*force.y;
            this->velocity.z += timestep*force.z;
        }

        Point velocity;
        Point position;

        int id;
        double mass;
};
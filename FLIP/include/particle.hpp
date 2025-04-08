
#define NUM_PARTICLES 1024
#define SQRT_NUM_PARTICLES 32

struct Point
{
    public:
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
};

struct Particle
{
    public:
        Particle(double x, double y, double z)
        {
            // TODO make sure this is memory safe in terms of ownership
            this->velocity = Point();
            this->position = Point(x, y, z);
        };

        Particle()
        {
            // TODO make sure this is memory safe
            this->velocity = Point();
            this->position = Point();
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

        double mass;
};
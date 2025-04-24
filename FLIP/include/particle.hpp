#include <sstream>
#include <string>
#include <algorithm>

#define NUM_PARTICLES 2000

/*

    Position coordinate system:

    y - image height
    ^
    |
    |                   image width
    |                   |
    .-----------------> x
    (0, 0)
*/

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

    // todo -> fix this to be bounded within cell size
    double ngp_distance(Point p1, double width, double height)
    {
        // printf("p1.x: %f, p1.y: %f\n", p1.x, p1.y);
        // printf("p2.x: %f, p2.y: %f\n", x, y);
        double dist = abs(x - p1.x)<=width/2 && abs(y - p1.y)<=height/2 ? 1 : 0;
        return dist;
    }

    // TODO -> this is not right (1 - (x - xi))
    double bilinear(Point p1, double width, double height, double depth)
    {
        double hdist = (1.0 - abs(this->x - p1.x)/width)*(1.0 - abs(this->y - p1.y)/height);
        double ddist = hdist * (1.0 - abs(this->z - p1.z)/depth);
        return std::clamp(ddist, 0.0, 1.0);
    }

    std::string print() {
        std::stringstream ss;
        ss << "p.x "<<this->x<<", p.y "<< this->y<<"\n";
        return ss.str();
    };
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

        std::string print() {return position.print();};

        Point velocity;
        Point position;

        int id;
        double mass;
};
#include <map>

#define PROFILE_MODE 1
#define DEBUG 0
#define REALTIME 0
#define VERBOSE 1

// how many timesteps to use for profiling
#define NUM_TIME_STEPS 100

// WARNING: THE MACROS ASSUME THE CLASS YOU ARE PROFILING HAS A POINTER TO THE PROFILER NAMED AS this->profiler

// only define profile_fun_start to do something if profile_mode is defined
#if PROFILE_MODE
#define PROFILE_FUN_START this->profiler->profile_start(__FUNCTION__)
#else
#define PROFILE_FUN_START
#endif

// only define profile_fun_end to do something if profile_mode is defined
#if PROFILE_MODE
#define PROFILE_FUN_END this->profiler->profile_end(__FUNCTION__)
#else
#define PROFILE_FUN_END
#endif

// only define profile_start to do something if profile_mode is defined
#if PROFILE_MODE
#define PROFILE_START(name) this->profiler->profile_start(name)
#else
#define PROFILE_START(_name)
#endif

// only define profile_end to do something if profile_mode is defined
#if PROFILE_MODE
#define PROFILE_END(name) this->profiler->profile_end(name)
#else
#define PROFILE_END(_name)
#endif

// only define profile_print to do something if profile_mode is defined
#if PROFILE_MODE
#define PROFILE_PRINT this->profiler->profile_print()
#else
#define PROFILE_PRINT
#endif

typedef std::chrono::time_point<std::chrono::steady_clock> time_start_t;
typedef std::chrono::duration<double> duration_t;

/*
A performance profiler for the particle simulator, used to analyze parts of the program runtime
*/
class Profiler{
    public:
        Profiler()
        {
            this->function_profiler = {};
            this->profiler_scratch = {};
        };

        void profile_start(std::string key);
        void profile_end(std::string key);
        void profile_print();

    private:
        std::map<std::string, duration_t> function_profiler;
        std::map<std::string, time_start_t> profiler_scratch;
};

void Profiler::profile_start(std::string name)
{
    this->profiler_scratch[name] = std::chrono::steady_clock::now();
}

void Profiler::profile_end(std::string name){
  duration_t delta = (std::chrono::steady_clock::now() - this->profiler_scratch[name]);
  if(this->function_profiler.count(name)>0){
    this->function_profiler[name] = this->function_profiler[name] + delta;
  }
  else this->function_profiler[name] = delta;
}

void Profiler::profile_print() {
  printf("Profiling results: \n");
  for (const auto& [key, value] : this->function_profiler){
    double compute_time = std::chrono::duration_cast<duration_t>(value).count();
    printf("[%s] = %f sec \n", key.c_str(), compute_time);
  }
}
#include <lang/Thread.h>

namespace lang {

::pthread_mutex_t Mutex::_mutexInitializer = PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;

}

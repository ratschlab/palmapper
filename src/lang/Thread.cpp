#include <stdio.h>
#include <lang/Thread.h>

namespace lang {

::pthread_mutex_t Mutex::_mutexInitializer = PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;

int Thread::launch() {
   return ::pthread_create(&_handle, NULL, &launchInterface, this);
}

void *Thread::launchInterface(void *self) {
   try {
      ((Thread*) self)->run();
   } catch (...) {
      ::printf("Unhandled exception in thread\n");
   }
   return NULL;
}

}

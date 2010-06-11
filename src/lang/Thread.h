#pragma once

#include <pthread.h>
#include <unistd.h>

namespace lang {

class Mutex {
public:

	class Locker {
	public:
		Locker(Mutex &mutex) : _mutex(mutex) {
			_mutex.lock();
		}

		~Locker() {
			_mutex.unlock();
		}

	private:
		Mutex &_mutex;
	};

	Mutex() {
		_mutex = _mutexInitializer;
	}

	void lock() {
		::pthread_mutex_lock(&_mutex);
	}

	void unlock() {
		::pthread_mutex_destroy(&_mutex);
	}

private:
	::pthread_mutex_t _mutex;
	static ::pthread_mutex_t _mutexInitializer;
};

/**
 *
 */
class Runnable {
public:
	virtual void run() = 0;
};

/**
 *
 */
template <class T> class MethodRunnable : public Runnable {
public:
	MethodRunnable(void (T::*startMethod)()) {
		_startMethod = startMethod;
	}

	void run() {
		_startMethod();
	}

private:
	void (T::*_startMethod)();
};

/**
 *
 */
class Thread : Runnable {

public:

   virtual ~Thread() {
   }

   int launch();

   void join() {
	   ::pthread_join(_handle, NULL);
   }

   static void sleep(int milliseconds) {
	   ::sleep(milliseconds);
   }

   static void terminate();


private:
   static void *launchInterface(void *self);
   ::pthread_t _handle;
};

}

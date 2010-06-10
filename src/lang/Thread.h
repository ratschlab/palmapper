#pragma once

#include <pthread.h>

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

class Thread {
public:

};

}

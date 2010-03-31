#pragma once

struct STORAGE_ENTRY{
	unsigned char id[4];
};

struct INDEX_ENTRY {
	unsigned int num;
	//STORAGE_ENTRY *last_entry;
	unsigned int offset;
};


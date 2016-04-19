// textfile.h
// --- interface for reading and writing text files
//
// simple reading and writing for text files
//
// www.lighthouse3d.com
//
// You may use these functions freely.
// they are provided as is, and no warranties, either implicit,
// or explicit are given
//////////////////////////////////////////////////////////////////////


#ifndef Text_File_h
#define Text_File_h

char *textFileRead(char *fn);
int textFileWrite(char *fn, char *s);

#endif
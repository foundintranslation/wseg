// wseg.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <cstdio>
#include <iostream>
#include <cassert>
#include <time.h>
#include "File.h"
#include "option.h"
#include "Vocab.h"
#include "Ngram.h"
#include "MLSeg.h"

static char *lmFile="small.1bo";
//"/n/gosling/transitory/mhwang/LM/expt1/BN.p1e-10.2bo.gz";
//"/g/ssli/data/mandarin/RT-04-Mandarin-CTS/scripts/TN/ML/trigram.txt";

static char *testData=NULL, *outfilename=NULL;
static unsigned int order = 1;
static unsigned debug=2, skipn=0;
int keep_skip=0, keepsym=0;
static unsigned maxchar=1000000, garbage_len=2000;
static int removespace=0;
static double unkpen=1e-10, prune=1e-50;
static int lcase=0, LFM=0, skipoov=0;
static int ucase=0;
static char *unk_in_history="@reject@";
static char *freqfile=0;
static char *symfile=NULL, *removefile=NULL;
static char *comments = "##";
static unsigned int nbest=1;

static Option options [] = {
	{OPT_TRUE, "keepsym", &keepsym, "keep punctuations and symbols (unused now)"},
	{OPT_STRING, "freq", &freqfile, "freq file"},
	{OPT_STRING, "lm", &lmFile, "file in ARPA LM format"},
	{OPT_STRING, "mapoov", &unk_in_history, "map OOV in history to this vocab"},
	{OPT_STRING, "in", &testData, "test text"},
	{OPT_STRING, "out", &outfilename, "output segmented text file"},
	{OPT_STRING, "comments", &comments, "comment lines begin with this"},
	{OPT_STRING, "symfile", &symfile, "symbol mapping file"},
	{OPT_STRING, "removefile", &removefile, "lines containing these characters will be removed"},
	{OPT_UINT, "order", &order, "ngram order"},
	{OPT_UINT, "nbest", &nbest, "nbest search"},
	{OPT_UINT, "debug", &debug, "lm debug level"},
	{OPT_UINT, "skipn", &skipn, "skip first n fields each line"},
	{OPT_TRUE, "keep_skip", &keep_skip, "keep the skipped fields in output"},
	{OPT_UINT, "maxchar", &maxchar, "max#chars per word"},
	{OPT_TRUE, "lcase", &lcase, "convert non-ch words into lower case"},
	{OPT_TRUE, "ucase", &ucase, "convert non-ch words into upper case"},
	{OPT_TRUE, "LFM", &LFM, "longest-first match"},
	{OPT_TRUE, "removespace", &removespace, "remove space between chs before seg"},
	{OPT_FLOAT, "unkPen", &unkpen, "unk penalty"},
	{OPT_FLOAT, "prune", &prune, "pruning threshold for nbest seg"},
	
	// These two parameters controled by main()
	{OPT_TRUE, "skipoov", &skipoov, "Skip sentences that have OOV words"},
	{OPT_UINT, "garbage", &garbage_len, "If input too long, throw it away as garbage"},

};

int main(int argc, char **argv)
{
	Vocab *vocab;
	Ngram *ngramLM;
	LM *useLM;
	Boolean ok;
	int argc1;
	char **argv1;

	if (0)
	{ 
		FILE *fp = fopen("tmp.txt", "w");
		unsigned char t[20];
		unsigned short code=0xa2a4, c;
		int i;
		
		t[2] = 0;
		c = code;
		for (i=0; i < 200; i++) {
			t[0] = (c & 0xFF00) >> 8;
			t[1] = c & 0xFF;
			fprintf(fp, "%s  space\n", t);
			c++;
		}
		fclose(fp);
		exit(0);
	}

	if (argc <= 1) {
		printf("option error\n");
		return 0;
	}

	if (argc == 2 && argv[1][0] == '@') {
     FILE *fp = fopen(&argv[1][1], "r");  
     int len;
     char word[256], *cptr;

     if (fp==NULL) {
	fprintf(stderr, "Can't open arg file %s\n", &argv[1][1]);
        exit(-1);
     }
     argc1 = 1;
     len=strlen(argv[0])+1;
     while (fscanf(fp, "%s ", word)==1) { len += strlen(word)+1; argc1++;}
     rewind(fp);
     cptr = new char[sizeof(char *) * argc1 + len]; 
     argv1 = (char **) cptr;
     cptr = (char *) (argv1 + argc1);

     argc1 = 0;
     argv1[argc1++] = cptr;
     strcpy(cptr, argv[0]);
     cptr += strlen(argv[0]) + 1;
     while (fscanf(fp, "%s ", word)==1) {
	argv1[argc1++] = cptr;
	strcpy(cptr, word);
	cptr += strlen(word)+1;
     }
     fclose(fp);
  } else {
     argc1 = argc;
     argv1 = argv;
  }

  Opt_Parse(argc1, argv1, options, Opt_Number(options), 0);
  if (argv != argv1) {
      delete [] argv1;
  }
  assert(lmFile);
  if (LFM) order = 1;

  vocab = new Vocab;
  assert(vocab);
  vocab->unkIsWord() = true; 
  /* <unk> is always added into vocab.
     But still may not be in LM.
  */
  {
       File file(lmFile, "r");

       ngramLM = new Ngram(*vocab, order);
       assert(ngramLM);
       ngramLM->debugme(debug);

       ok = ngramLM->read(file, false);
       assert(ok);
       useLM = ngramLM;
  }

  if (LFM) 
	printf("Longest first match from lexicon %s\n",
           lmFile);
  else
     printf("numWords=%d, %d-gram will be used from %s\n",
        vocab->numWords(), order, lmFile);
  fflush(stdout);

  if (!testData) return 0;

  int transform = NO_TOUCH;
  time_t s1, s2;
  clock_t t1, t2;
  FILE *file;
  FILE *outfp=stdout;
  unsigned numsent=0, numchars=0;
  LogP total_logp=0;
  int total_segs=0, nseg, len, lno=0;
  double p0, p1;
  Result out;
  char *line;
  cMLSeg Segment(vocab, useLM, order, debug, comments);

  if ((line = new char[garbage_len+1]) == NULL)
    exit(1);

  if (ucase) transform = UPPER;
  if (lcase) transform = LOWER;
  if (ucase && lcase) {
         cerr << "-lcase and -ucase are exclusive\n";
         return (-1);
  }

  if ((file=fopen(testData, "r")) == NULL) {
           cerr << "Can't open " << testData << endl;
           exit(-1);
  }

  if (outfilename) {
    outfp = fopen(outfilename, "w");
    if (outfp == NULL) {
           cerr << "Can't create " << outfilename << endl;
           return -1;
    }
  }

  Segment.SetParam(maxchar, transform, removespace, unkpen, LFM);
  Segment.SetParam2(skipn, keep_skip);
  Segment.SetParam3(nbest, prune);
  Segment.SetUnkInHis(unk_in_history);
  Segment.LoadMap(symfile);
  if (Segment.LoadRemoved(removefile) < 0) exit(-1);
  if (Segment.InitializeMem() < 0) exit(-1);
    
  s1=time(NULL); 
  t1=clock();

  while (fgets(line, garbage_len+1, file)) {
    lno++;
    len = strlen(line);
    assert(len <= garbage_len);
    /* sometimes the input file does not end with \n. 
       But the last line is not garbage.
    */
    if (line[len-1] != '\n' && len == garbage_len) {
      int c;
      fprintf(stderr, "%s(%d) seems to be garbage (%s)\n", 
              testData, lno, line);
      while ((c=fgetc(file)) != '\n' && c != EOF);

      fprintf(outfp, "\n");
      continue;
    }
    if (line[len-1] == '\n') 
	line[len-1] = 0; /* so that a garbage gb byte at the end won't be combined with\n */
    else {
	fprintf(stderr, "last line does not end with newline (%s)\n", line);
    }

    nseg = Segment.segment1(line, &out);

    if (nseg <= 0) {
      // maybe a comment line, or a line with a garbage GB code
      if (out.str) fprintf(outfp, "%s\n", out.str);
      else fprintf(outfp, "\n");

      if (debug >= 4)
        fprintf(stderr, "%s(%d) nseg=%d\n", testData, lno, nseg);

      Segment.release1();
      continue;
    }


    assert(out.str);
    if (skipoov && out.oov) 
      fprintf(outfp, "\n");
    else  {
      total_logp += out.score; // includes score to </s>
      total_segs += nseg; // not include </s>
      numchars += out.chars;
      numsent++;

      fprintf(outfp, "%s\n", out.str);
      if (debug >= 3 && !LFM) {
        p0 = LogPtoPPL(out.score/(nseg+1));
        p1 = LogPtoPPL(out.score/nseg);
        printf("\tlogprob=%g, ppl=%g, ppl1=%g\n", out.score, p0, p1);
      }

      if ((numsent & 0x0FFF) == 0)  // every 4k utterances
         cout << numsent << " valid sentences processed." << endl << endl;

      if ((numsent & 0x0FFFFF) == 0) // every 1M utterances
        Segment.DumpStatistics();
    }

    Segment.release1();
  }

  fclose(file);
  if (outfp != stdout) fclose(outfp);

  if (freqfile)
    Segment.PrintOOV(freqfile, true);
  else if (debug >= 2)
    Segment.PrintOOV(NULL, false);

  if (numsent) {
      printf("\nAmong non-empty output sentences:\n");
      printf("num_sent=%d, num_words = %d, num_chars=%d\n",
             numsent, total_segs, numchars);

      if (!LFM) {
        p0 = LogPtoPPL(total_logp/(total_segs+numsent));
        p1 = LogPtoPPL(total_logp/total_segs);

        printf("Word ppl=%g, ppl1=%g logprob=%g\n",
               p0, p1, total_logp);
    
        p0 = LogPtoPPL(total_logp/(numchars+numsent));
        p1 = LogPtoPPL(total_logp/numchars);
        printf("Char ppl=%g, ppl1=%g\n", p0, p1);
      }
   }

   t2=clock();
   s2=time(NULL);
   //float sec=(t2-t1)/CLOCKS_PER_SEC,h;
   float sec=(s2-s1), h;
   h=sec/3600;
   printf("\n%.3f secs (%.3f hours) taken for segmenting %d lines\n", 
         sec, h, lno);

   delete line;

   return 0;
}

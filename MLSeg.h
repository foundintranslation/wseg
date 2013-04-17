#pragma once
#ifndef _SEARCH_H
#define _SEARCH_H
#include "Vocab.h"
#include "LM.h"

extern const char *GBmap[];

#define MAX_STATE   500000
#define MAX_HISTORY 500000

#define SPACE	" @sp@ "

#define NO_TOUCH 0
#define UPPER    1
#define LOWER    2

struct stat_ml {
   unsigned max;
   char *maxcase;
};


struct Result {
   char *str;
   LogP score;
   int nseg, chars, oov;
};

struct History {
   VocabIndex wid;
   LogP score;
   History *backptr;
   unsigned short fromi, len;
};

struct State {
   History *history;
   union {
     struct {
       unsigned short w0, w1;
     };
     unsigned int biword_his;
   };

   union {
     struct {
       unsigned short fromi, len;
     };
     unsigned int inputpos;
   };
};

class cMLSeg {
  public:
  char *m_savedstr;
     void SetParam2(unsigned skip, bool keepskip) {
       m_skip = skip;
       m_keepskip = keepskip;
     }

     void SetParam3(unsigned nbest, double prune) {
       m_nbest = nbest;
       m_prune = log10(prune);
     }

     void SetParam(unsigned maxchar, int transform, bool removespace, 
            double unkPen, bool LFM=false)
     {
        m_maxchar = maxchar;
        m_transform = transform;
        m_keepspace = !removespace;
        m_unkpen = log10(unkPen);
        m_LFM = LFM;
     }


     // Pr(w|...OOV...) = Pr(w|...mapunk...)
     void SetUnkInHis (char *mapunk) {
       if (!mapunk) return;

       m_oov_in_his = m_vocab->getIndex(mapunk);
       if (m_oov_in_his == Vocab_None) {
         cerr << mapunk << " not in vocab" << endl;
         exit(-1);
       }
     }


     cMLSeg(Vocab *v, LM *m, unsigned order, unsigned debug, char *comment) 
     {
        m_comment = comment;
        m_GBmap = (char **)GBmap;
        m_Removed = NULL;
        m_skip = 0;
        m_vocab = v;
        m_LM = m;
        m_orderm1 = order - 1;
        m_transform = NO_TOUCH;
        m_keepspace = true;
        m_verbose = debug;
        m_vocsize = v->numWords();

        m_unk = m_oov_in_his = m_vocab->unkIndex();

        m_LFM = false;
        m_freq = NULL;
        m_fsize = 0;

        m_totalstates = m_totalframes = 0;
        m_totalstates_i = m_totalframes_i = 0;

        memset(&m_maxhis, 0, sizeof(stat_ml));
        memset(&m_maxstates, 0, sizeof(stat_ml));
        memset(&m_maxout, 0, sizeof(stat_ml));
        m_unkpen = log10(1e-10); // Pr(OOV|h)= m_unkpen
        m_prune = log10(1e-50); // pruning threshold for nbest search
        m_nbest = 1;
        m_outstring = m_outend = m_savedstr = NULL;
        m_inIndex = NULL;
        m_CurrentState = m_NextState = NULL;
        m_History = m_border = NULL;
        m_Secmap = NULL;
     }


     ~cMLSeg() {
        DumpStatistics();

        delete [] m_CurrentState;
        delete [] m_NextState;
        delete [] m_freq;
        delete [] m_History;

        if (m_Secmap)
          delete m_Secmap[0];
        delete m_Secmap;

        delete m_Removed;

        free(m_maxhis.maxcase);
        free(m_maxstates.maxcase);
        free(m_maxout.maxcase);
     }


     int InitializeMem();
     char **LoadMap(char *filename);
     int LoadRemoved(char *filename);

     int segment1(VocabString instring, Result *outstring);
     unsigned PrintOOV(char *filename, bool all);

     void DumpStatistics() {
           cout << endl;

           if (m_maxstates.max) 
             printf("maxstates at one frame = %d for (%s)\n", m_maxstates.max,
                m_maxstates.maxcase);

           if (m_totalframes || m_totalframes_i) {
             double states = m_totalstates + m_totalstates_i;
             double frames = m_totalframes + m_totalframes_i;
                       
             printf("avg #states/frame = %.0f/%.0f = %.2f\n",
                    states, frames, states/frames);
           }

           if (m_maxhis.max) 
             printf("maxhis=%d for (%s)\n", m_maxhis.max, m_maxhis.maxcase);

           if (m_maxout.max) 
	      printf("longest output (%d bytes) = %s\n",
	         m_maxout.max, m_maxout.maxcase);
     }

     void release1() {
          delete [] m_outstring;
          delete [] m_savedstr;
          delete [] m_inIndex;
     };
 
  private:
     bool IsOOV(VocabIndex v) {
          return v==Vocab_None ||
                  v==m_unk || (v >= m_vocsize);
     };

     int ReplaceChars (const char *in, char **outp);
     char *Find1 (char **map, const char *iptr, char *optr);
     char *Find2 (char **map, const char *iptr, char *optr);
     int normalize_text(char *src);
     char *Add1Entry(char **map, char *src, char *dst, char *bptr);
     bool BeRemoved(char *sptr);

     int CheckCode(char *in);
     int MLSearch(unsigned starti, unsigned endi);
     int LFMSearch(unsigned starti, unsigned endi);
     void InitializeSearch1(unsigned starti);
     int BacktraceAnswer(unsigned endi);
     unsigned AddNextState(State *s);
     VocabIndex GetPhraseIndex(unsigned start, unsigned nw);
     void DumpInput(unsigned int start);
     void DumpStates(State *s, unsigned len, bool dumphis);
     unsigned IncrementFreq(VocabIndex i);

     void UpdateStatistics(struct stat_ml *s, unsigned n) {
         if (n > s->max) {
            s->max = n;
            free(s->maxcase);
            s->maxcase = strdup(m_savedstr);
         }
     }  

     Vocab *m_vocab;
     LM *m_LM;
     VocabIndex m_unk, m_oov_in_his;
     /* m_unk = m_vocab->unkIndex(). Just a quick way to get

        Pr(OOV|h) = m_unkpen = -big s.t. OOV is discouraged.

        Pr(w|...OOV...) = Pr(w|...oov_in_his...)
        where m_oov_in_his is usually @reject@.
     */ 
     unsigned int m_vocsize, m_verbose, m_time;
     unsigned *m_freq, m_fsize;

     double m_totalstates, m_totalframes;
     int m_totalstates_i, m_totalframes_i;
     /* Whenever int is overflown, put the big number in float and
        let the int keeps the left-over.
        So totalstates = m_totalstates + m_totalstates_i;
     */

     char **m_Secmap, **m_GBmap, *m_Removed;
     LogP m_unkpen, m_prune;
     struct stat_ml m_maxhis, m_maxstates, m_maxout;
     unsigned  m_maxchar, m_orderm1, m_nbest;
     char *m_comment;
     bool m_keepspace, m_keepskip, m_LFM;
     int m_skip, m_transform;

     // For each input sentence
     LogP m_finalscore;
     History *m_finalpath;
     char *m_outstring, *m_outend, **m_inIndex;
     unsigned m_in, m_chars, m_oov;

     // search space
     unsigned m_hindex, m_cindex, m_nindex;
     State *m_CurrentState, *m_NextState;
     History *m_History, *m_border;
};

#endif

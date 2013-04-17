#include "stdafx.h"
#include <iostream>
#include "File.h"
#include <assert.h>
#include <ctype.h>
#include <locale>
#include "option.h"
#include "Vocab.h"
#include "Ngram.h"
#include "MLSeg.h"

int strincmp (const char *str1, const char *str2, int n);
int nChineseChar(unsigned char *s);


int cMLSeg::InitializeMem()
{
  if (m_orderm1 != 0 && m_orderm1 != 1 && m_orderm1 != 2) {
      cerr << "Current works for only unigram, bigram, or trigram ML seg!" 
	<< endl;

      return -2;
  }

  if (m_vocsize >= 65536) {
      cerr << "History prunning won't work! VocSize too big!" << endl;
      return -2;
  }

  m_History = new History[MAX_HISTORY];
  m_CurrentState = new State[MAX_STATE];
  m_NextState = new State[MAX_STATE];
  m_freq = new unsigned int [m_vocsize];
  
  if (!m_CurrentState || !m_NextState || !m_freq || !m_History) {
      fprintf(stderr, "%s(%d): running out of memory\n",
        __FILE__, __LINE__);
      return -1;
  }

  m_border = &m_History[2]; 
  /* m_Histroy[0] not used. Reserved for future purposes.
     m_History[1] = <s>
     Always starts real history from m_History[2] = m_border
  */
  memset(m_History, 0, sizeof(History)* (m_border-m_History));
  memset(m_freq, 0, sizeof(unsigned) * m_vocsize);
  m_fsize = m_vocsize;

  VocabIter viter(*m_vocab, false);
  VocabString s, maxs;
  int longestchar=0, err=0;
  int chsw=0;

  while ((s=viter.next())) 
      if (*s &0x80) 
      {
          int len = strlen(s);
          int nchar = nChineseChar((unsigned char *)s);

          if (len != nchar*2) {
               cerr << "mixed English+Chinese? " << s << endl;
               err++;
          } else {
              if (nchar > m_maxchar) {
                cerr << s << " in vocab longer than " << m_maxchar << " chars" 
	         << endl;
              } else {
                chsw++;
                if (nchar > longestchar) {
	          longestchar = nchar;
                  maxs = s;
                }
              }
           }
      }

    /* If user allows unlimited long chs word (m_maxchar=infinite),
       but our longestchar is only n (small), then update m_maxchar.
    */
    if (m_maxchar > longestchar)
       m_maxchar = longestchar;

  if (longestchar) {
    cout << chsw << " Chinese words will be used in segmentation." << endl;
    printf("longest word = %s (%d)\n", maxs, m_maxchar);
  } else {
    return -3;
  }

  if (err) return -err;

  return 0;
}



void cMLSeg::InitializeSearch1(unsigned starti)
{
  unsigned i;
  History *his;
  VocabIndex w0;
  
  his = m_border-1;
  memset(his, 0, sizeof(History));
  his->wid = m_vocab->ssIndex(); // <s>

  m_CurrentState->history = his;
  m_CurrentState->fromi = starti;
  m_CurrentState->len = 0;
  m_CurrentState->biword_his = 0xFFFFFFFF;
  if (m_orderm1) // bigram or trigram
    m_CurrentState->w1 = his->wid; // <s>

  m_cindex = 1;
  m_hindex = m_border - m_History;
  m_finalpath = NULL;
  m_finalscore = 0;
}



int cMLSeg::segment1(VocabString in, Result *out)
{
   char *mapin=NULL;
   unsigned start;
   int total_nseg;

   m_savedstr = m_outstring = m_outend = NULL;
   m_inIndex = NULL;
   m_chars = m_oov = m_in = 0;

   if (!in || !out || *in==0) {
	total_nseg = -2;
        goto cleanup;
   }

   if (strncmp(in, m_comment, strlen(m_comment)) == 0) {
     int len=strlen(in);
     m_outstring = new char[len+1];
     m_outend = m_outstring + len;
     assert(m_outstring); 
     strcpy(m_outstring, in);

     if (m_outstring[len-1] == '\n')
	m_outstring[len-1] = 0;
     total_nseg = 0;

     goto cleanup;
   }

   ReplaceChars(in, &mapin);
   if ((total_nseg=CheckCode(mapin)) < 0 )
      goto cleanup;

   start = 0;
   assert(m_inIndex[start][0] != ' ');  // maybe 0

   if (m_skip) {
      char *src, c;
      unsigned  f, i;

      for (i=f=0; f < m_skip && i < m_in; i++) 
         if (m_inIndex[i][0] == ' ') {
             f++;
         } 

      start = i;
      if (m_keepskip) {
        src = m_inIndex[i];
        c = *src;
        *src = 0;
        strcpy(m_outstring, m_inIndex[0]);
        *src = c;
      }
   }
   assert(m_inIndex[start][0] != ' '); // maybe 0

   if (m_verbose==100) {
     DumpInput(start);
     total_nseg = 0;
     goto cleanup;
   }

   if (start >= m_in) { 
      total_nseg = 0;
      goto cleanup;
   }
  
   if (m_LFM)
     total_nseg = LFMSearch(start, m_in);
   else {
     InitializeSearch1(start);
     total_nseg = MLSearch(start, m_in);
   }

cleanup:
   assert(m_outstring==NULL || m_outstring+strlen(m_outstring) <= m_outend);
   delete [] mapin;

   out->str = m_outstring;
   out->score = m_finalscore;
   out->oov = m_oov;
   out->chars = m_chars;
   out->nseg = total_nseg;

   return out->nseg;
}



void cMLSeg::DumpInput(unsigned int start)
{
  unsigned i;
  char *aptr, *bptr, c;

  printf("%s\n", m_savedstr);
  for (i=0; i < m_in; i++) {
    aptr = m_inIndex[i];
    bptr = m_inIndex[i+1];
    c = *bptr;
    *bptr = 0;
    if (i == start) printf("*");
    printf("(%d): %s\n", i, aptr);
    *bptr = c;
  }
  printf("\n");
}


unsigned cMLSeg::PrintOOV(char *filename, bool all)
{
   unsigned int i, start, countoov=0;
   unsigned int limit=m_vocab->numWords(); // including OOVs
   VocabString s;
   FILE *fp=stdout;

   assert(m_fsize >= limit);

   cout << "\n";
   if (filename) {
     if ((fp = fopen(filename, "w")) == NULL) {
       fprintf(stderr, "Can't create %s\n", filename);
       return -1;
     }
   } 

   start = all? 0 : m_vocsize;

   for (i=start; i < limit; i++) 
     if (m_freq[i]) {
       s = m_vocab->getWord(i);
       fprintf(fp, "%-20s %5d ", s, m_freq[i]);
       
       if (i >= m_vocsize) {
          countoov++;
          fprintf(fp, "OOV-%d ", countoov);
       }
       fprintf(fp, "\n");
     }

   if (fp != stdout) fclose(fp);

   cout << countoov << " OOV words" << endl;

   return countoov;
}


unsigned cMLSeg::IncrementFreq(VocabIndex i)
{
  if (i >= m_fsize) {
    unsigned newsize = m_fsize + i + 1000;
    unsigned *newf = new unsigned[newsize];
    unsigned j;

    assert(newf);
    memcpy(newf, m_freq, sizeof(unsigned) * m_fsize);
    for (j=m_fsize; j < newsize; j++)
      newf[j] = 0;

    delete [] m_freq;
    m_freq = newf;
    m_fsize = newsize;
  }

  m_freq[i]++;
  if (IsOOV(i)) m_oov++;

  return 1;
}


VocabIndex cMLSeg::GetPhraseIndex(unsigned start, unsigned nw)
{
  char *eptr, c;
  VocabIndex index;

  if (start+nw > m_in) 
     return Vocab_None;

  eptr = m_inIndex[start+nw];
  c = *eptr;
  *eptr = 0;
  index = m_vocab->getIndex(m_inIndex[start]);
  *eptr = c;

  return index;   // maybe Vocab_None
}



/* Return negative if error. 0 otherwise
   O/P: m_inIndex[0] never begins with a white space, but
        may end with a white space.
*/
int cMLSeg::CheckCode(char *in)
{
  char *dptr, *cptr;
  unsigned len = strlen(in), end_with_space;
  locale loc;
  /* From in to m_savedstr:
     worst case: Each byte adds two white spaces as an English word does.
  */
  len = len*3;
  m_savedstr = new char[len+1]; // for 0 at the end
  m_inIndex = new char *[len+1];
  m_outstring = new char [2*len+1];
  /* From m_savedstr to m_outstring.
     In the worst case, every byte in m_savedstr is inserted a white space
     after it.
  */

  assert(m_savedstr);
  assert(m_inIndex);
  assert(m_outstring);
  *m_outstring = 0;
  m_outend = m_outstring + 2*len;

  cptr = in;
  dptr = m_savedstr;
  m_in = m_chars = 0;
  while (isspace(*cptr,loc)) cptr++;

  if (m_skip) {
    for (int i=0; i < m_skip && *cptr; i++) {
      m_inIndex[m_in++] = dptr;

      for (; !isspace(*cptr,loc) && *cptr; cptr++) 
        *dptr++ = *cptr;

      m_inIndex[m_in++] = dptr;
      *dptr++ = ' ';

      while (isspace(*cptr,loc)) cptr++;
    }
  }

  assert(! isspace(*cptr,loc));
  end_with_space = 1; // so that m_inIndex[] won't begin with a space

  while (*cptr) {
    if (*cptr & 0x80) {
        m_inIndex[m_in++] = dptr;
        if (BeRemoved(cptr))
          return -4;

        *dptr++ = *cptr++;
        *dptr++ = *cptr;
        if (*cptr == 0) {
           fprintf(stderr, "missing one more byte at EOL\n");
           return -3;
        }

        cptr++;
        m_chars++;
        end_with_space = 0;
        continue;
     }

     // *cptr = white space
     if (isspace(*cptr, loc)) {
       do { cptr++; } while (isspace(*cptr, loc));

       if (m_keepspace && !end_with_space) {
           m_inIndex[m_in++] = dptr;
           *dptr++ = ' ';

           end_with_space = 1;
       }
       continue;
     }

     // *cptr = English word: insert a white space before and after
     if (! end_with_space) {
       m_inIndex[m_in++] = dptr;
       *dptr++ = ' ';

       end_with_space = 1;
     }

     m_inIndex[m_in] = dptr;
     do {
       *dptr++ = *cptr++;
     } while (*cptr && (*cptr&0x80)==0 && !isspace(*cptr,loc));
     *dptr = 0;

     if (strcmp(m_inIndex[m_in], "@sp@")==0) {
       /* don't need to keep @sp@ in the input.
          " @sp@ " will force white space to replace @sp@
          because it is an English word.
       */
       dptr = m_inIndex[m_in];
       continue;
     }

     m_chars++; // 1 English word counted as 1 char
     if (m_transform != NO_TOUCH && strcmp(m_inIndex[m_in], "[laugh]") != 0 &&
         strcmp(m_inIndex[m_in], "@reject@") != 0)
     // upper/lower case unification
     { char *xp;
       for (xp=m_inIndex[m_in]; xp < dptr; xp++)
         if (*xp >= 'a' && *xp <= 'z') {
           if (m_transform == UPPER) *xp += 'A' - 'a';
         } else if (*xp >= 'A' && *xp <= 'Z') {
           if (m_transform == LOWER) *xp += 'a' - 'A';
         }
     }
     end_with_space = 0;
     m_in++;

     m_inIndex[m_in++] = dptr;
     *dptr++ = ' ';
     end_with_space = 1;

     continue;
  }

  /*
  if (end_with_space && m_in) {
    dptr = m_inIndex[--m_in];
    assert(*dptr == ' ');
  }
  */

  *dptr = 0;
  m_inIndex[m_in] = dptr; // trick for GetPhraseIndex() 

  assert(dptr <= m_savedstr + len);
  assert(m_in <= len);
 
  return 0;
}



unsigned cMLSeg::AddNextState(State *s)
{
  State *a;
  unsigned i;

   // Viterbi pruning for speedup.
   // If top n segmentations are asked, can't do this.
   if (m_nbest==1) {
     for (i=0,a=m_NextState; i < m_nindex; i++, a++) 
       {
         if (a->biword_his == s->biword_his &&
             a->inputpos == s->inputpos)
           {
             if (s->history->score > a->history->score) 
               *a = *s;
             return 0;
           }
       }
   } else {
     for (i=0,a=m_NextState; i < m_nindex; i++, a++) 
       {
         if (a->biword_his == s->biword_his &&
             a->inputpos == s->inputpos)
           {
             if (s->history->score < (a->history->score + m_prune))
               return 0;
           }
       }
   }

   if (m_nindex >= MAX_STATE) {
	fprintf(stderr, "NextState out of mem (t=%d/%d) for (%s)\n", 
                m_time, m_in, m_savedstr);
        exit(-1);
   }
   m_NextState[m_nindex++] = *s;
   return 1;
}



void cMLSeg::DumpStates(State *statearray, unsigned slen, bool dumphis)
{
   unsigned int i;
   State *s;
   History *h;
   char *ptr, *bptr, c;

   cout << "t=" << m_time << endl;
   s = statearray;
   for (i=0; i < slen; i++, s++) {
     bptr = m_inIndex[s->fromi];
     ptr=m_inIndex[s->fromi+s->len];
     c = *ptr;
     *ptr = 0;
     printf("s=%d input=%s, h%d\n", i, bptr, s->history - m_History);
     *ptr = c;
   }

   if (dumphis) {
     h = m_History;
     for (i=0; i < m_hindex; i++, h++) {
         fprintf(stdout, "h%d (%s %.2f h%d)\n", i,
	    h->wid==0xFFFF ? "null" : m_vocab->getWord(h->wid), 
            h->score,
            h->backptr ? (h->backptr - m_History) : -1);
     }
   }
   cout << endl;
}

int cMLSeg::LFMSearch(unsigned starti, unsigned endi)
{
   char *optr, *bptr, *eptr, c;
   int beginning, remaining_len, len, index, nw;

   optr = m_outstring + strlen(m_outstring);
   beginning = starti;
   remaining_len = endi - starti;

   nw = 0;
   while (remaining_len) {   
     bptr = m_inIndex[beginning];
     if (*bptr == ' ') {
       beginning++;
       remaining_len--;
       continue;
     }

     len = remaining_len;
     index = GetPhraseIndex(beginning, len);
     while (IsOOV(index) && len != 1) {
       len--;
       index = GetPhraseIndex(beginning, len);
     }

     eptr = m_inIndex[beginning + len];
     c = *eptr;
     *eptr = 0;

     strcpy(optr, bptr);
     optr += eptr - bptr;
     *optr++ = ' ';

     index = m_vocab->addWord(bptr);
     IncrementFreq(index);
     *eptr = c;

     nw++;
     beginning += len;
     remaining_len -= len;
   }
   *optr = 0;

   m_finalscore = nw; 
   return nw;
}



int cMLSeg::MLSearch(unsigned starti, unsigned endi)
{
   unsigned i, j, y, lastchannel;
   State *state, newstate, *ns;
   VocabIndex index, lmindex, history[20];
   History *his;
   LogP logp;
   int sum;

   sum = m_totalframes_i + (endi - starti);
   if (sum <= 0) {
     m_totalframes += m_totalframes_i;
     m_totalframes_i = endi - starti;
   } else
     m_totalframes_i = sum;

   if (m_verbose >= 4) {
	char *bptr = m_inIndex[starti];
        char *eptr = m_inIndex[endi];
        char c = *eptr;
        *eptr = 0;
	printf("%s\n", bptr);
        *eptr = c;
   }

   for (m_time=starti; m_time < endi; m_time++) {
      m_nindex = 0;
      if (m_verbose >= 4)
         DumpStates(m_CurrentState, m_cindex, true);

      if (m_inIndex[m_time][0] == ' ') {
         // force a word boundary at the previous char
         state = m_CurrentState;
         for (i=j=0; i < m_cindex; i++, state++)
           if (state->len == 0) {
               m_CurrentState[j] = *state;
               m_CurrentState[j++].fromi++;
           }
         m_cindex = j;

         // count #states AFTER each frame m_time.
         UpdateStatistics (&m_maxstates, m_cindex);
         sum = m_totalstates_i + m_cindex;
         if (sum <= 0) {
           m_totalstates += m_totalstates_i;
           m_totalstates_i = m_cindex;
         } else 
           m_totalstates_i = sum;

         continue;
      }

      // Add one more char to all partial words
      for (j=0; j < m_cindex; j++) {
        state = &m_CurrentState[j];
        if (state->len >= m_maxchar) 
           continue;
        
        state->len++;
        m_NextState[m_nindex++] = *state;
      }

      // Now check word boundary
      lastchannel = m_nindex;
      for (j=0; j < lastchannel; j++) {
         state = &m_NextState[j];
         index = GetPhraseIndex(state->fromi, state->len);
         lmindex = index;

         if (!IsOOV(index) || state->len==1) 
         {
            if (IsOOV(index)) {
	      // OOV is guaranteed to be single-char or an English word
              // Either way, state->len==1
                assert(state->len==1);
                lmindex = m_oov_in_his;
                logp = m_unkpen;
            } else {
                his = state->history;
                for (y=0; his && y < m_orderm1; y++) {
                  history[y] = his->wid;
                  his = his->backptr;
                }
                history[y] = Vocab_None;
                logp = m_LM->wordProb(lmindex, history);
            }

            if (m_hindex >= MAX_HISTORY) {
              fprintf(stderr, "History out of mem (t=%d/%d) for (%s)\n", 
                      m_time, endi, m_savedstr);
               exit(-1);
            }
            his = &m_History[m_hindex++];
            his->wid = lmindex;
            his->fromi = state->fromi;
            his->len = state->len;
            his->score = state->history->score + logp;
	    his->backptr = state->history;

            ns = &newstate;
            ns->history = his;
            ns->fromi = state->fromi + state->len;
            assert (ns->fromi == m_time+1);
            ns->len = 0;
            ns->biword_his = 0xFFFFFFFF; // unigram
            if (m_orderm1 >= 1) // 2- or 3-gram
              ns->w1 = lmindex;
            if (m_orderm1 == 2) // 3-gram
              ns->w0 = his->backptr->wid;
            
            AddNextState(ns);
         } 
      } // for j

      UpdateStatistics (&m_maxstates, m_nindex);
      sum = m_totalstates_i + m_nindex;
      if (sum <= 0) { // overflow
        m_totalstates += m_totalstates_i;
        m_totalstates_i = m_nindex;
      } else 
        m_totalstates_i = sum;

      ns = m_CurrentState;
      m_CurrentState = m_NextState;
      m_NextState = ns;
      m_cindex = m_nindex;
   }

   UpdateStatistics(&m_maxhis, m_hindex);

   return BacktraceAnswer(endi);
}



int cMLSeg::BacktraceAnswer(unsigned endi)
{
   State *bs, *state;
   unsigned i;
   int y, nw;
   char *cptr, *eptr, *optr, c;
   VocabIndex *wid, history[20];
   History *his;
   LogP logp;
   bool sent_end = (endi >= m_in);

   if (m_verbose >= 4)
      DumpStates(m_CurrentState, m_cindex, true);

   bs = NULL;
   for (i=0; i < m_cindex; i++) {
      state = &m_CurrentState[i];
      if (state->len != 0) continue;

      logp = state->history->score;

      if (sent_end) {
        // now transit to </s>
        VocabIndex seindex = m_vocab->seIndex();
        his = state->history;
        for (y=0; his && y < m_orderm1; y++) {
           history[y] = his->wid;
           his = his->backptr;
        }
        history[y] = Vocab_None;
        logp += m_LM->wordProb(seindex, history);
      } 

      if (!bs || logp > m_finalscore)
      {
          bs = state;
          m_finalscore = logp;
      }
   }

   if (bs==NULL) {
      fprintf(stderr, "bug?");
      exit(-1);
   }
   assert(bs);
   m_finalpath = bs->history;

   for (his=bs->history, y=0; his >= m_border; his=his->backptr) 
	y++;
   assert(y);

   wid = new VocabIndex[y];
   assert(wid);
   nw = y;
   for (his=bs->history; his >= m_border; his=his->backptr)
   {
       // his->wid can be unk. So not usable here. Has to use string
       cptr = m_inIndex[his->fromi];
       eptr = m_inIndex[his->fromi + his->len];
       c = *eptr;
       *eptr = 0;
       //printf("%s %.5f\n", cptr, his->score - his->backptr->score);

       wid[--y] = m_vocab->addWord(cptr); // add OOV if necessary
       *eptr = c;
   }
   assert(y==0);

   optr = m_outstring + strlen(m_outstring);
   for (y=0; y < nw; y++) {
        IncrementFreq(wid[y]);
	cptr = (char *) m_vocab->getWord(wid[y]);
        strcpy(optr, cptr);
	optr += strlen(cptr);
	*optr++ = ' ';
   }
   *optr = 0;
   delete [] wid;

   UpdateStatistics(&m_maxout, strlen(m_outstring));

   return nw; // does not include </s>
}


char *cMLSeg::Find2 (char **map, const char *iptr, char *optr) {
   int i;
   char *src, *dst;

   for (i=0; map[i]; i +=2) {
      src = map[i];
      if (*src == *iptr && *(src+1) == *(iptr+1)) {
	 dst = map[i+1];
	 strcpy(optr, dst);
	 optr += strlen(dst);
         return optr;
      }
    }

    return NULL;
}


char *cMLSeg::Find1 (char **map, const char *iptr, char *optr) {
   int i;
   char *src, *dst;

   for (i=0; map[i]; i += 2) {
      src = map[i];
      if (*src == *iptr) {
	 dst = map[i+1];
	 strcpy(optr, dst);
         optr += strlen(dst);
         return optr;
      }
    }

    return NULL;
}

            
int cMLSeg::ReplaceChars (const char *in, char **outp) {
   const char *iptr;
   char *optr, *op;
   int count=0, len;

   len = strlen(SPACE)*strlen(in) + 1;
   *outp = new char[len];
   assert(*outp);
   optr = *outp;

   for (iptr=in; *iptr; iptr++) {
 	if (*iptr & 0x80) {
           op = Find2(m_Secmap, iptr, optr);
	   if (op==NULL) 
	      op = Find2(m_GBmap, iptr, optr);
	   
	   if (op) {
	     optr = op;
	     count++;
           } else {
              *optr++ = *iptr;
	      *optr++ = *(iptr+1);
	   }

	   iptr++;
	   continue;
        } 
	
        op = Find1(m_Secmap, iptr, optr);
	if (op==NULL)
 	    op = Find1(m_GBmap, iptr, optr);

	if (op) {
	  optr = op;
	  count++;
	} else {
	    *optr++ = *iptr;
        }
   }
   *optr = 0;

   assert(optr < *outp + len);
// printf("in=%s, out=%s (%d)\n", in, *outp, count);

   return count;
}



char **cMLSeg::LoadMap (char *filename) {
  FILE *fp;
  int n, bytes, m;
  char line[256], src[80], dst[80], *bptr;

  if (filename==NULL || (fp=fopen(filename, "r")) == NULL) {
    if (filename) fprintf(stderr, "Can't open %s\n", filename);

    n = 2;
    m_Secmap = new char*[2*n+1];
    bytes = (3+strlen(SPACE)+1) * 2;
    m_Secmap[0] = NULL;
    bptr = new char[bytes];

    goto addws;
  }

  n = bytes = 0;
  while (fgets(line, sizeof(line)-1, fp) != NULL) {
    if (sscanf(line, "%s %s", src, dst) != 2) continue;

    bytes += strlen(src)+1;
    if (strcmp(dst, "space") == 0) bytes += strlen(SPACE)+1;
    else bytes += strlen(dst)+1;
    n++;
  }
  printf("%d entries in %s\n", n, filename);
  n += 2; // add 0xa88c and 0x9c9b
  bytes += (3 + strlen(SPACE)+1) * 2;

  m_Secmap = new char* [2*n+1];
  bptr = new char[bytes];
  assert(m_Secmap);
  m_Secmap[0] = NULL;

  rewind(fp);
  while (fgets(line, sizeof(line)-1, fp) != NULL) {
    if (sscanf(line, "%s %s", src, dst) != 2) continue;

    bptr = Add1Entry(m_Secmap, src, dst, bptr);
    if (bptr == NULL) exit(-1);
  }
  fclose(fp);

addws:
  strcpy(src, "0xa88c");
  strcpy(dst, "space");
  bptr = Add1Entry(m_Secmap, src, dst, bptr);

  strcpy(src, "0x9c9b");
  strcpy(dst, "space");
  bptr = Add1Entry(m_Secmap, src, dst, bptr);
  assert(bptr <= m_Secmap[0] + bytes);

  for (m=0; m_Secmap[m]; m++);
  assert(m <= 2*n);

  return m_Secmap;
}



int cMLSeg::normalize_text(char *src)
{
  if (strcmp(src, "space") == 0) {
    strcpy(src, SPACE);
    return 2;
  }

  if (strncmp(src, "0x", 2) == 0 && strlen(src)==6) {
    int i;
    char a, h, l, save[20];

    strcpy(save, src);
    for (i=2; i < 6; i++) {
      a = src[i];
      if (a >= 'A' && a <= 'Z') a += 'a' - 'A';

      if (a >= 'a' && a <= 'z') 
        src[i] = (a - 'a') + '9' + 1;
    }

    h = (src[2]-'0') * 16 + (src[3]-'0');
    l = (src[4]-'0') * 16 + (src[5]-'0');
    src[0] = h;
    src[1] = l;
    src[2] = 0;

    return 1;
  }

  return 0;
}
    


char *cMLSeg::Add1Entry(char **map, char *src, char *dst, char *bptr)
{
  int i;

  normalize_text(src);
  normalize_text(dst);

  for (i=0; map[i]; i += 2)
    if (strcmp(map[i], src) == 0) 
    {
        if (strcmp(map[i+1], dst) == 0)
          return bptr;
        else {
          fprintf(stderr, "diff mapping for %s => %s vs. %s\n", 
               src, map[i+1], dst);
          return NULL;
        }  
    }

  map[i++] = bptr;
  strcpy(bptr, src);
  bptr += strlen(src)+1;

  map[i++] = bptr;
  strcpy(bptr, dst);
  bptr += strlen(bptr)+1;

  map[i] = NULL;
  return bptr;
}


bool cMLSeg::BeRemoved(char *srcptr) {

  for (char *bptr=m_Removed; *bptr; bptr += 2) {
    if (*bptr==*srcptr && *(bptr+1)==*(srcptr+1)) 
      return true;
  }

  return false;
}

int cMLSeg::LoadRemoved(char *filename)
{
  FILE *file;
  int n;
  char line[256], word[256], *bptr;

  if (filename == NULL) {
    m_Removed = new char[1];
    m_Removed[0] = 0;
    return 0;
  }

  if ((file=fopen(filename, "r")) == NULL) {
    fprintf(stderr, "Can't open %s\n", filename);
    return -1;
  }

  n = 0;
  while (fgets(line, sizeof(line)-1, file)) 
    if (sscanf(line, "%s", word) == 1) n++;
  n++;
  m_Removed = new char[2*n+1];
  assert(m_Removed);

  rewind(file);
  bptr = m_Removed;
  n = 0;
  while (fgets(line, sizeof(line)-1, file)) 
    if (sscanf(line, "%s", word) == 1 &&
        strlen(word)==2) 
    {
      if (word[0] & 0x80) {
        *bptr++ = word[0];
        *bptr++ = word[1];
        n++;
      }
    }
  fclose(file);
  *bptr = 0;
  assert(bptr <= m_Removed + 2*n);

  printf("%d removed chars in %s\n%s\n", n, filename, m_Removed);
  return n;
}



int strincmp (const char *str1, const char *str2, int n)
{
   int i, a, b;

   for (i=0; i < n; i++) {
      a = str1[i];
      b = str2[i];
      if (a >= 'a' && a <= 'z') a += 'A' - 'a';
      if (b >= 'a' && b <= 'z') b += 'A' - 'a';
      if (a!=b) return a-b;
      if (a==0 && b==0) return 0;
   }
   return 0;
}

int nChineseChar(unsigned char *s)
{
   int c = 0;
   unsigned char *t;

   for (t=s; *t; t++) {
      if (*t >= 128) {
        c++;
        t++;
      }
   }

   return c;
}


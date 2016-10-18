/*
 nw_cmdline.c
 project: NeedlemanWunsch
 author: Isaac Turner <turner.isaac@gmail.com>
 url: http://sourceforge.net/projects/needlemanwunsch
 Copyright (C) 06-Dec-2011
 
 see: README

 == License
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h> // tolower
#include <stdarg.h> // required for va_list

// my utility functions
#include "string_buffer.h"
#include "bioinf.h"
#include "utility_lib.h"

// Alignment scoring and loading
#include "alignment_scoring_load.h"

#include "needleman_wunsch.h"

// For this run
char* cmd;
char print_colour = 0, print_pretty = 0, print_scores = 0,
     print_fasta = 0, print_zam = 0;

SCORING_SYSTEM* scoring = NULL;

// Alignment results stored here
char *alignment_a = NULL, *alignment_b = NULL;
t_buf_pos alignment_max_length;


void set_default_scoring()
{
  scoring = scoring_system_default();
}

void print_usage(char* err_fmt, ...)
{
  if(err_fmt != NULL)
  {
    va_list argptr;
    va_start(argptr, err_fmt);
    
    STRING_BUFFER *error = string_buff_init(200);
    string_buff_append_str(error, "NeedlemanWunsch Error: ");
    string_buff_vsprintf(error, string_buff_strlen(error), err_fmt, argptr);
    string_buff_chomp(error);

    va_end(argptr);

    fprintf(stderr, "%s\n", error->buff);
    fprintf(stderr, "Use -h option to print help\n");
    exit(EXIT_FAILURE);
  }

  if(scoring != NULL)
  {
    scoring_free(scoring);
  }

  // Get and print defaults
  set_default_scoring();

  fprintf(stderr, "usage: %s [OPTIONS] [seq1 seq2]\n", cmd);

  fprintf(stderr,
"  Needleman-Wunsch optimal global alignment (maximises score). Takes a pair \n"
"  of sequences on the command line, reads from a file and from sequence \n"
"  piped in.  Can read gzip files and those in FASTA, FASTQ or plain format.\n"
"\n"
"  OPTIONS:\n"
"    --file <file>        Sequence file reading with gzip support\n"
"    --files <f1> <f2>    Read one sequence from each file at a time to align\n"
"    --stdin              Read from STDIN (same as '--file -')\n"
"    --case_sensitive     Case sensitive character comparison\n"
"\n"
"    --scoring <PAM30|PAM70|BLOSUM80|BLOSUM62>\n"
"    --substitution_matrix <file>  see details for formatting\n"
"    --substitution_pairs <file>   see details for formatting\n"
"\n");

  fprintf(stderr,
"    --match <score>      [default: %i]\n"
"    --mismatch <score>   [default: %i]\n"
"    --gapopen <score>    [default: %i]\n"
"    --gapextend <score>  [default: %i]\n",
          scoring->match, scoring->mismatch,
          scoring->gap_open, scoring->gap_extend);

  fprintf(stderr,
"\n"
"    --freestartgap       No penalty for gap at start of alignment\n"
"    --freeendgap         No penalty for gap at end of alignment\n"
"\n"
"    --printscores        Print optimal alignment scores\n"
"    --printfasta         Print fasta header lines\n"
"    --pretty             Print with a descriptor line\n"
"    --colour             Print with colour\n"
"    --zam                A funky type of output\n"
"\n"
" DETAILS:\n"
"  * For help choosing scoring, see the README file. \n"
"  * Gap (of length N) penalty is: (open+N*extend)\n"
"  * To do alignment without affine gap, set '--gapopen 0'.\n"
"  * Scoring files should be matrices, with entries separated by a single \n"
"    character or whitespace.  See files in the 'scores' directory for examples.\n"
"\n"
"  turner.isaac@gmail.com  (compiled: "COMPILE_TIME")\n");

  exit(EXIT_FAILURE);
}

void align_zam(char *seq_a, char *seq_b)
{
  needleman_wunsch(seq_a, seq_b, alignment_a, alignment_b, scoring);

  // Swap '-' for '_'
  int i;
  for(i = 0; alignment_a[i] != '\0'; i++)
  {
    if(alignment_a[i] == '-')
    {
      alignment_a[i] = '_';
    }

    if(alignment_b[i] == '-')
    {
      alignment_b[i] = '_';
    }
  }

  int num_of_mismatches = 0;
  int num_of_indels = 0;

  // Print branch 1
  printf("Br1:%s\n", alignment_a);

  // Print spacer
  printf("    ");

  for(i = 0; alignment_a[i] != '\0'; i++)
  {
    if(alignment_a[i] == '_' || alignment_b[i] == '_')
    {
      printf(" ");
      num_of_indels++;
    }
    else if((scoring->case_sensitive && alignment_a[i] != alignment_b[i]) ||
            tolower(alignment_a[i]) != tolower(alignment_b[i]))
    {
      printf("*");
      num_of_mismatches++;
    }
    else
    {
      printf("|");
    }
  }

  printf("\n");

  // Print branch 2
  printf("Br2:%s\n", alignment_b);

  // print mismatch indel numbers
  printf("%i %i\n\n", num_of_mismatches, num_of_indels);
}

void align(char *seq_a, char *seq_b,
           char *seq_a_name, char *seq_b_name)
{
  if(print_zam)
  {
    return align_zam(seq_a, seq_b);
  }

  int score = needleman_wunsch(seq_a, seq_b, alignment_a, alignment_b, scoring);

  if(print_fasta && seq_a_name != NULL)
  {
    printf("%s\n", seq_a_name);
  }

  if(print_fasta && print_pretty && seq_b_name != NULL)
  {
    printf("%s\n", seq_b_name);
  }

  if(print_colour)
  {
    // Print alignment line 1
    alignment_colour_print_against(alignment_a, alignment_b,
                                   scoring->case_sensitive);
  }
  else
  {
    printf("%s", alignment_a);
  }
  printf("\n");
  
  if(print_pretty)
  {
    // Print spacer
    alignment_print_spacer(alignment_a, alignment_b, scoring);

    printf("\n");
  }
  else if(print_fasta && seq_b_name != NULL)
  {
    printf("%s\n", seq_b_name);
  }

  if(print_colour)
  {
    // Print alignment line 2
    alignment_colour_print_against(alignment_b, alignment_a,
                                   scoring->case_sensitive);
  }
  else
  {
    printf("%s", alignment_b);
  }
  printf("\n");

  if(print_scores)
  {
    printf("score: %i\n", score);
  }
  
  printf("\n");
}

// If seq2 is NULL, read pair of entries from first file
// Otherwise read an entry from each
void align_from_file(SEQ_FILE *seq1, SEQ_FILE *seq2)
{
  STRING_BUFFER *entry1_title = string_buff_init(200);
  STRING_BUFFER *entry2_title = string_buff_init(200);
  STRING_BUFFER *entry1_seq = string_buff_init(200);
  STRING_BUFFER *entry2_seq = string_buff_init(200);

  char empty_file = 1;

  while(1)
  {
    seq_file_read(seq1, entry1_title, entry1_seq);

    if(string_buff_strlen(entry1_seq) == 0)
    {
      break;
    }
    else
    {
      empty_file = 0;
    }

    seq_file_read((seq2 == NULL ? seq1 : seq2), entry2_title, entry2_seq);

    if(string_buff_strlen(entry2_seq) == 0)
    {
      fprintf(stderr, "Odd number of sequences - I read in pairs!\n");
      break;
    }

    // Align
    char *title1 = NULL, *title2 = NULL;

    if(seq_file_get_type(seq1) != SEQ_PLAIN)
    {
      title1 = entry1_title->buff;
    }

    if((seq2 == NULL && seq_file_get_type(seq1) != SEQ_PLAIN) ||
       (seq2 != NULL && seq_file_get_type(seq2) != SEQ_PLAIN))
    {
      title2 = entry2_title->buff;
    }

    // Check memory
    t_buf_pos new_max_alignment = entry1_seq->len + entry2_seq->len;

    if(new_max_alignment > alignment_max_length)
    {
      // Expand memory used for storing result
      alignment_max_length = new_max_alignment;

      if(!nw_realloc_mem((unsigned int)new_max_alignment,
                         &alignment_a, &alignment_b))
      {
        print_usage("Ran out of memory");
      }
    }

    align(entry1_seq->buff, entry2_seq->buff, title1, title2);
  }

  if(empty_file)
  {
    fprintf(stderr, "NeedlemanWunsch: Warning, empty input\n");
  }

  // Free memory
  string_buff_free(entry1_title);
  string_buff_free(entry2_title);
  string_buff_free(entry1_seq);
  string_buff_free(entry2_seq);
}


double align1(char* seq_a, char* seq_b, int * c, int * m)
{
  // Variables to store alignment result
  char *alignment_a, *alignment_b;

  // malloc the above variables
  // (seq1 and seq2 are used to figure out how much memory may be needed)
  nw_alloc_mem(seq_a, seq_b, &alignment_a, &alignment_b);

  // Decide on scoring
  int match = 1;
  int mismatch = -2;
  int gap_open = -4;
  int gap_extend = -1;
  
  // Don't penalise gaps at the start
  // ACGATTT
  // ----TTT would score +3 (when match=+1)
  char no_start_gap_penalty = 0;
  
  // ..or gaps at the end e.g.
  // ACGATTT
  // ACGA--- would score +4 (when match=+1)
  char no_end_gap_penalty = 0;

  // Compare character case-sensitively (usually set to 0 for DNA etc)
  char case_sensitive = 0;

  SCORING_SYSTEM* scoring = scoring_create(match, mismatch,
                                           gap_open, gap_extend,
                                           no_start_gap_penalty,
                                           no_end_gap_penalty,
                                           case_sensitive);

  // Add some special cases
  // x -> y means x in seq1 changing to y in seq2
  scoring_add_mutation(scoring, 'a', 'c', -2); // a -> c give substitution score -2
  scoring_add_mutation(scoring, 'c', 'a', -1); // c -> a give substitution score -1

  // We could also prohibit the aligning of characters not given as special cases
  // scoring->use_match_mismatch = 0;

  needleman_wunsch(seq_a, seq_b, alignment_a, alignment_b, scoring);
  int matchCount1 = 0;
  int nCount = 0;
  int i = 0;
  while(i < strlen(alignment_b) && i < strlen(alignment_a))
  {
      if(alignment_a[i] == alignment_b[i])
      {
          matchCount1++;
      }
      if(alignment_b[i] == 'N')
      {
          nCount++;
      }else if(alignment_b[i] != '-' && alignment_a[i] != alignment_b[i])
      {
          (*m)++;
      }
      
      i++;
  }
  
  // Free memory used to store scoring preferences
  scoring_free(scoring);
  free(alignment_a);
  free(alignment_b);
  *c = *c + matchCount1;
  double rr = (double)(matchCount1)/(double)(strlen(seq_a));
  return rr;
  
}


char * GetContigFromContigSet(char * contigSetFile, long int contigIndex){
    
    char * resultContig = NULL;
    
    long int maxSize = 10000;
    char * contig = NULL;
    if(NULL == (contig = (char*)malloc(sizeof(char)*maxSize))){
        perror("malloc error!");
        exit(1);
    }
    
    FILE * fp; 
    if((fp = fopen(contigSetFile, "r")) == NULL){
        printf("%s, does not exist!", contigSetFile);
        exit(0);
    }
    
    long int contigCount = -1;
    
    while((fgets(contig, maxSize, fp)) != NULL){ 
       
       if(contig[0] == '>'){  
           contigCount++;
           if(contigCount > contigIndex){
               break;
           }
           continue;
           
       }
       if(contigCount == contigIndex){
           
           long int extendLength = strlen(contig);
           if(contig[extendLength-1] == '\n'){
               extendLength--;
           }
           long int contigLength = 0;
           char * tempContig = NULL;
           if(resultContig != NULL){
               contigLength = strlen(resultContig);
               tempContig = (char *)malloc(sizeof(char)*(contigLength+1));
               strncpy(tempContig, resultContig, contigLength);
               free(resultContig);
                   
               resultContig = (char *)malloc(sizeof(char)*(contigLength + extendLength + 1));
               strncpy(resultContig, tempContig, contigLength);
                       
               strncpy(resultContig + contigLength, contig, extendLength);
               resultContig[contigLength + extendLength] = '\0';
               free(tempContig);
           }else{
               resultContig = (char *)malloc(sizeof(char)*(extendLength+1));
               strncpy(resultContig, contig, extendLength);
               resultContig[extendLength] = '\0';
           }       
       }
    }  
    
    fclose(fp);
    
    long int i = 0;
    while(i < strlen(resultContig)){
        resultContig[i] = toupper(resultContig[i]);
        i++;
    }
    
    return resultContig;
}

int main(int argc, char* argv[])
{
  
  int count = 0;
  int i = 0;
  int count9 = 0;
  int count8 = 0;
  int count5 = 0;
  int a = 0;
  int b = 0;
  int * c = &a;
  int * m = &b;
  long int gapLength = 0;
  while(i < atoi(argv[3])){
      
      char * left = GetContigFromContigSet(argv[1], i);
      char * right = GetContigFromContigSet(argv[2], i);
      gapLength = gapLength + strlen(left);
      double r = align1(left, right, c, m);
      
      if(r >= atof(argv[4])){
          count++;
      }
      if(r >= 0.9){
          count9++;
      }
      if(r >= 0.8){
          count8++;
      }
      if(r >= 0.5){
          count5++;
      }
      
      i++;
  }
  
  FILE * fp;
  fp = fopen(argv[5], "a+");
  fprintf(fp, "gap_count:%i\n", i);
  fprintf(fp, "all_gap_length:%ld\n", gapLength);
  fprintf(fp, "match_count:%i\n", a);
  fprintf(fp, "mismatch_count:%i\n", b);

  double p = (double)a/(double)(a+b);
  double r = (double)a/(double)gapLength;
  double f = 2*p*r/(p+r);
  fprintf(fp, "precision:%f\n", p);
  fprintf(fp, "recall:%f\n", r);
  fprintf(fp, "f-score:%f\n\n", f);
  
  return 1;
}

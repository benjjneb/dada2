#include <stdlib.h>
#include "dada.h"
// #if _WIN32
// #include <intrin.h>
// #endif
// [[Rcpp::interfaces(cpp)]]

/************* KMERS *****************
 * Current Kmer implementation assumes A/C/G/T only.
 * Current Kmer implementation requires kmer size of 3-7 only.
 */

double kmer_dist(uint16_t *kv1, int len1, uint16_t *kv2, int len2, int k) {
  int i;
  int n_kmer = 1 << (2*k); // 4^k kmers
  uint16_t dotsum = 0;
  double dot = 0.0;
  
  for(i=0;i<n_kmer; i++) {
    dotsum += (kv1[i] < kv2[i] ? kv1[i] : kv2[i]);
  }
  
  dot = ((double) dotsum)/((len1 < len2 ? len1 : len2) - k + 1.);
  
  return (1. - dot);
}

#ifdef __x86_64 //  Computes kmer distance with SSE2 intrinsics.
double kmer_dist_SSEi(uint16_t *kv1, int len1, uint16_t *kv2, int len2, int k) {
  size_t n_kmer = 1 << (2*k); // 4^k kmers
  int16_t dst[64];
  size_t STEP=8;
  uint16_t dotsum = 0;
  double dot = 0.0;
  
  __m128i vsum = _mm_set1_epi16(0);
  
  for(uint16_t const * end( kv1 + n_kmer );kv1<end;kv1+=STEP,kv2+=STEP) { // assumes n_kmer divisible by STEP
    __m128i kk1 = _mm_loadu_si128( (__m128i*) kv1 );
    __m128i kk2 = _mm_loadu_si128( (__m128i*) kv2 );
    __m128i kmin = _mm_min_epi16(kk1, kk2);
    vsum = _mm_add_epi16(vsum, kmin);
  }
  _mm_storeu_si128 ((__m128i*) &dst, vsum);
  dotsum = dst[0] + dst[1] + dst[2] + dst[3] + dst[4] + dst[5] + dst[6] + dst[7];
  
  dot = ((double) dotsum)/((len1 < len2 ? len1 : len2) - k + 1.);
  return (1. - dot);
}
#else // not X64, define dummy functions that won't be called
double kmer_dist_SSEi(uint16_t *kv1, int len1, uint16_t *kv2, int len2, int k) {
  return 0.0;
}
#endif

#ifdef __x86_64 //  Computes kmer distance with SSE2 intrinsics.
// Uses kmers packed into 8 bit integers. Can overflow, in which case negative value is returned.
double kmer_dist_SSEi_8(uint8_t *kv1, int len1, uint8_t *kv2, int len2, int k) {
  size_t n_kmer = 1 << (2*k); // 4^k kmers
/// memcpy to put the two vectors in spatial proximity solves the cache miss issue and speed everything up greatly.
/// However it then becomes the bottleneck.
///  uint8_t kv[2048];
///  memcpy(kv, kv1, 1024);
///  memcpy(&kv[1024],kv2,1024);
///  kv1=kv;
///  kv2=&kv[1024];
  uint8_t dst[64];
  size_t STEP=16;
  int i;
  uint16_t dotsum = 0;
  double dot = 0.0;
  
  __m128i vsum_a = _mm_set1_epi8(0x00);
//  __m128i vsum_b = _mm_set1_epi8(0x00);
  
  for(uint8_t const * end( kv1 + n_kmer );kv1<end;kv1+=STEP,kv2+=STEP) { // assumes n_kmer divisible by STEP
    __m128i kk1 = _mm_loadu_si128( (__m128i*) kv1 );
    __m128i kk2 = _mm_loadu_si128( (__m128i*) kv2 );
    __m128i kmin_a = _mm_min_epu8(kk1, kk2);
    vsum_a = _mm_adds_epu8(vsum_a, kmin_a);
  }
  _mm_storeu_si128 ((__m128i*) &dst, vsum_a);

  bool overflow=false;
  for(i=0;i<STEP;i++) {
    if(dst[i] == 255) { overflow=true; }
    dotsum += dst[i];
  }
  if(overflow) { return(-1.); }
  
  dot = ((double) dotsum)/((len1 < len2 ? len1 : len2) - k + 1.);
  return (1. - dot);
}
#else // not X64, define dummy functions that won't be called
double kmer_dist_SSEi_8(uint8_t *kv1, int len1, uint8_t *kv2, int len2, int k) {
  return 0.0;
}
#endif

// Computes kmer distance based on ordered kmers (e.g. gapless align)
// If different lengths, returns -1 (invalid)
double kord_dist(uint16_t *kord1, int len1, uint16_t *kord2, int len2, int k) {
  int i;
  uint16_t dotsum = 0;
  double dot = 0.0;
  
  if(len1 != len2 || kord1 == NULL || kord2 == NULL) { return(-1.0); } // Exit if different lengths
  size_t klen = len1 - k + 1;

  for(i=0;i<klen;i++) {
    dotsum += (kord1[i] == kord2[i]);
  }
  
  dot = ((double) dotsum)/((len1 < len2 ? len1 : len2) - k + 1.);
  return (1. - dot);
}

#ifdef __x86_64
// Computes kmer distance with SSE2 intrinsics, based on ordered kmers (e.g. gapless align)
// If different lengths, returns -1 (invalid)
double kord_dist_SSEi(uint16_t *kord1, int len1, uint16_t *kord2, int len2, int k) {
  int i;
  int16_t dst[64];
  size_t STEP=8;
  uint16_t dotsum = 0;
  double dot = 0.0;
  
  if(kord1 == NULL || kord2 == NULL) { return(-1.0); }
  size_t klen = (len1 < len2 ? len1 : len2) - k + 1;
  size_t n_vec = klen - (klen % STEP); // Vectorize over this many
  
  __m128i vsum = _mm_set1_epi16(0);
  
  for(uint16_t const * end( kord1 + n_vec );kord1<end;kord1+=STEP,kord2+=STEP) {
    __m128i kk1 = _mm_loadu_si128( (__m128i*) kord1 );
    __m128i kk2 = _mm_loadu_si128( (__m128i*) kord2 );
    __m128i keq = _mm_cmpeq_epi16(kk1, kk2);
    vsum = _mm_sub_epi16(vsum, keq); // Using 0xFFFF = -1
  }
  _mm_storeu_si128 ((__m128i*) &dst, vsum);
  for(i=0;i<STEP;i++) {
    dotsum += dst[i]; 
  }
  for(i=n_vec;i<klen;i++,kord1++,kord2++) { // kord starts pointing to where it was left
    if(*kord1 == *kord2) { dotsum++; };
  }
  
  dot = ((double) dotsum)/((len1 < len2 ? len1 : len2) - k + 1.);
  return (1. - dot);
}
#else // not X64, define dummy functions that won't be called
double kord_dist_SSEi(uint16_t *kord1, int len1, uint16_t *kord2, int len2, int k)
{
  return 0.0;
}
#endif

void assign_kmer8(uint8_t *kvec8, const char *seq, int k) {  // Assumes a clean seq (just 1s,2s,3s,4s)
  int i, j, nti;
  size_t len = strlen(seq);
  if(len <= 0 || len > SEQLEN) { Rcpp::stop("Unexpected sequence length."); }
  if(k >= len || k < 3 || k > 8) { Rcpp::stop("Invalid kmer-size."); }
  size_t klen = len - k + 1; // The number of kmers in this sequence
  size_t kmer = 0;
  size_t n_kmers = (1 << (2*k));  // 4^k kmers
  uint16_t *kvec = (uint16_t *) malloc(n_kmers * sizeof(uint16_t)); //E
  if (kvec == NULL)  Rcpp::stop("Memory allocation failed.");
  for(kmer=0;kmer<n_kmers;kmer++) { kvec[kmer] = 0; }
  
  if(len <=0 || len > SEQLEN) {
    Rcpp::stop("Unexpected sequence length.");
  }
  
  for(i=0; i<klen; i++) {
    kmer = 0;
    for(j=i; j<i+k; j++) {
      nti = ((int) seq[j]) - 1; // Change 1s, 2s, 3s, 4s, to 0/1/2/3
      if(nti != 0 && nti != 1 && nti != 2 && nti != 3) {
        Rcpp::stop("Unexpected nucleotide.");
        kmer = 999999;
        break;
      }
      kmer = 4*kmer + nti;
    }
    
    // Make sure kmer index is valid. This doesn't solve the N's/-'s
    // issue though, as the "length" of the string (# of kmers) needs
    // to also reflect the reduction from the N's/-'s
    if(kmer == 999999) { ; } 
    else if(kmer >= n_kmers) {
      Rcpp::stop("Kmer index out of range.");
    } else { // Valid kmer
      kvec[kmer]++;
    }
  }
  // Move to destination uint8_t vector w/ overflow checks
  for(kmer=0;kmer<n_kmers;kmer++) { 
    if(kvec[kmer]<255) {
      kvec8[kmer] = (uint8_t)kvec[kmer];
    } else {
      kvec8[kmer] = 255;
    }
  }
  free(kvec);
}

void assign_kmer(uint16_t *kvec, const char *seq, int k) {  // Assumes a clean seq (just 1s,2s,3s,4s)
  int i, j, nti;
  size_t len = strlen(seq);
  if(len <= 0 || len > SEQLEN) { Rcpp::stop("Unexpected sequence length."); }
  if(k >= len || k < 3 || k > 8) { Rcpp::stop("Invalid kmer-size."); }
  size_t klen = len - k + 1; // The number of kmers in this sequence
  size_t kmer = 0;
  size_t n_kmers = (1 << (2*k));  // 4^k kmers
  for(kmer=0;kmer<n_kmers;kmer++) { kvec[kmer] = 0; }
  
  if(len <=0 || len > SEQLEN) {
    Rcpp::stop("Unexpected sequence length.");
  }
  
  for(i=0; i<klen; i++) {
    kmer = 0;
    for(j=i; j<i+k; j++) {
      nti = ((int) seq[j]) - 1; // Change 1s, 2s, 3s, 4s, to 0/1/2/3
      if(nti != 0 && nti != 1 && nti != 2 && nti != 3) {
        Rcpp::stop("Unexpected nucleotide.");
        kmer = 999999;
        break;
      }
      kmer = 4*kmer + nti;
    }
    
    // Make sure kmer index is valid. This doesn't solve the N's/-'s
    // issue though, as the "length" of the string (# of kmers) needs
    // to also reflect the reduction from the N's/-'s
    if(kmer == 999999) { ; } 
    else if(kmer >= n_kmers) {
      Rcpp::stop("Kmer index out of range.");
    } else { // Valid kmer
      kvec[kmer]++;
    }
  }
}

// Returns the vector of the kmers in order
void assign_kmer_order(uint16_t *kord, char *seq, int k) {  // Assumes a clean seq (just 1s,2s,3s,4s)
  int i, j, nti;
  size_t len = strlen(seq);
  if(len <=0 || len > SEQLEN) { Rcpp::stop("Unexpected sequence length."); }
  if(k >= len || k < 1 || k > 8) { Rcpp::stop("Invalid kmer-size."); }
  size_t klen = len - k + 1; // The number of kmers in this sequence
  size_t kmer = 0;
  size_t n_kmers = (1 << (2*k));  // 4^k kmers
  if (kord == NULL)  Rcpp::stop("Memory allocation failed.");
  for(i=0;i<klen;i++) { kord[i] = 0; }
  
  for(i=0; i<klen; i++) {
    kmer = 0;
    for(j=i; j<i+k; j++) {
      nti = ((int) seq[j]) - 1; // Change 1s, 2s, 3s, 4s, to 0/1/2/3
      if(nti != 0 && nti != 1 && nti != 2 && nti != 3) {
        Rcpp::stop("Unexpected nucleotide.");
        kmer = 999999;
        break;
      }
      kmer = 4*kmer + nti;
    }
    
    // Make sure kmer index is valid. This doesn't solve the N's/-'s
    // issue though, as the "length" of the string (# of kmers) needs
    // to also reflect the reduction from the N's/-'s
    if(kmer == 999999) { ; } 
    else if(kmer >= n_kmers) { // maybe check if above uin16_t max? Defined in header <stdint.h>, UINT16_MAX
      Rcpp::stop("Kmer index out of range.");
    } else { // Valid kmer
      kord[i] = kmer;
    }
  }
}


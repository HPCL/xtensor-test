#include <iostream>
#include <time.h>
#include <memory>
#define ALIGN 16

constexpr size_t NN = 16*600000;//10000000;


void test_plainArray_matrix(float* const* input , const int type) {
  //
  float *Ax = new float[NN*36+4] __attribute__((aligned(ALIGN)));
  float *Bx = new float[NN*36+4] __attribute__((aligned(ALIGN)));
  float *Cx = new float[NN*36+4] __attribute__((aligned(ALIGN)));
  for (size_t j=0;j<NN*36;++j) Cx[j]=0.;
  //
  // store in matrix order (i.e. all elements of a given matrix are contiguous)
  for (size_t i=0;i<36;++i) {
    for (size_t x=0;x<NN;++x) {
      Ax[i + 36*x] = input[i][x];
      Bx[i + 36*x] = input[i][x]+0.5;
    }
  }

  const clock_t begin = clock();
#pragma vector aligned
  for (size_t x = 0; x < NN; ++x) {
    const size_t Nx = x*36;
    for (size_t i = 0; i < 6; ++i) {
      for (size_t j = 0; j < 6; ++j) {
        for (size_t k = 0; k < 6; ++k) {
          Cx[ Nx + (i*6 + j) ] += Ax[ Nx + (i*6 + k) ] * Bx[ Nx + (k*6 + j) ];
        }
      }
    }
  }
  const clock_t end = clock();

  // std::cout << "Ax=" << std::endl
  //        << Ax[(0*6+0)] << "\t" << Ax[(0*6+1)] << "\t" << Ax[(0*6+2)] << "\t" << Ax[(0*6+3)] << "\t" << Ax[(0*6+4)] << "\t" << Ax[(0*6+5)] << std::endl
  //        << Ax[(1*6+0)] << "\t" << Ax[(1*6+1)] << "\t" << Ax[(1*6+2)] << "\t" << Ax[(1*6+3)] << "\t" << Ax[(1*6+4)] << "\t" << Ax[(1*6+5)] << std::endl
  //        << Ax[(2*6+0)] << "\t" << Ax[(2*6+1)] << "\t" << Ax[(2*6+2)] << "\t" << Ax[(2*6+3)] << "\t" << Ax[(2*6+4)] << "\t" << Ax[(2*6+5)] << std::endl
  //        << Ax[(3*6+0)] << "\t" << Ax[(3*6+1)] << "\t" << Ax[(3*6+2)] << "\t" << Ax[(3*6+3)] << "\t" << Ax[(3*6+4)] << "\t" << Ax[(3*6+5)] << std::endl
  //        << Ax[(4*6+0)] << "\t" << Ax[(4*6+1)] << "\t" << Ax[(4*6+2)] << "\t" << Ax[(4*6+3)] << "\t" << Ax[(4*6+4)] << "\t" << Ax[(4*6+5)] << std::endl
  //        << Ax[(5*6+0)] << "\t" << Ax[(5*6+1)] << "\t" << Ax[(5*6+2)] << "\t" << Ax[(5*6+3)] << "\t" << Ax[(5*6+4)] << "\t" << Ax[(5*6+5)] << std::endl;
  // std::cout << "Bx=" << std::endl
  //        << Bx[(0*6+0)] << "\t" << Bx[(0*6+1)] << "\t" << Bx[(0*6+2)] << "\t" << Bx[(0*6+3)] << "\t" << Bx[(0*6+4)] << "\t" << Bx[(0*6+5)] << std::endl
  //        << Bx[(1*6+0)] << "\t" << Bx[(1*6+1)] << "\t" << Bx[(1*6+2)] << "\t" << Bx[(1*6+3)] << "\t" << Bx[(1*6+4)] << "\t" << Bx[(1*6+5)] << std::endl
  //        << Bx[(2*6+0)] << "\t" << Bx[(2*6+1)] << "\t" << Bx[(2*6+2)] << "\t" << Bx[(2*6+3)] << "\t" << Bx[(2*6+4)] << "\t" << Bx[(2*6+5)] << std::endl
  //        << Bx[(3*6+0)] << "\t" << Bx[(3*6+1)] << "\t" << Bx[(3*6+2)] << "\t" << Bx[(3*6+3)] << "\t" << Bx[(3*6+4)] << "\t" << Bx[(3*6+5)] << std::endl
  //        << Bx[(4*6+0)] << "\t" << Bx[(4*6+1)] << "\t" << Bx[(4*6+2)] << "\t" << Bx[(4*6+3)] << "\t" << Bx[(4*6+4)] << "\t" << Bx[(4*6+5)] << std::endl
  //        << Bx[(5*6+0)] << "\t" << Bx[(5*6+1)] << "\t" << Bx[(5*6+2)] << "\t" << Bx[(5*6+3)] << "\t" << Bx[(5*6+4)] << "\t" << Bx[(5*6+5)] << std::endl;
  std::cout << "Cx=" << std::endl
            << Cx[(0*6+0)] << "\t" << Cx[(0*6+1)] << "\t" << Cx[(0*6+2)] << "\t" << Cx[(0*6+3)] << "\t" << Cx[(0*6+4)] << "\t" << Cx[(0*6+5)] << std::endl
            << Cx[(1*6+0)] << "\t" << Cx[(1*6+1)] << "\t" << Cx[(1*6+2)] << "\t" << Cx[(1*6+3)] << "\t" << Cx[(1*6+4)] << "\t" << Cx[(1*6+5)] << std::endl
            << Cx[(2*6+0)] << "\t" << Cx[(2*6+1)] << "\t" << Cx[(2*6+2)] << "\t" << Cx[(2*6+3)] << "\t" << Cx[(2*6+4)] << "\t" << Cx[(2*6+5)] << std::endl
            << Cx[(3*6+0)] << "\t" << Cx[(3*6+1)] << "\t" << Cx[(3*6+2)] << "\t" << Cx[(3*6+3)] << "\t" << Cx[(3*6+4)] << "\t" << Cx[(3*6+5)] << std::endl
            << Cx[(4*6+0)] << "\t" << Cx[(4*6+1)] << "\t" << Cx[(4*6+2)] << "\t" << Cx[(4*6+3)] << "\t" << Cx[(4*6+4)] << "\t" << Cx[(4*6+5)] << std::endl
            << Cx[(5*6+0)] << "\t" << Cx[(5*6+1)] << "\t" << Cx[(5*6+2)] << "\t" << Cx[(5*6+3)] << "\t" << Cx[(5*6+4)] << "\t" << Cx[(5*6+5)] << std::endl;
  float time = float(end-begin)/CLOCKS_PER_SEC;
  std::cout << "plainArray_matrix -- time for NN=" << NN << " multiplications is " << time << " s, i.e. per track [s]=" << time/float(NN) << std::endl;

  delete Ax, Bx, Cx;
}

void test_plainArray_element(float* const* input, const int type) {
  float* Ax = new float[NN*36+4]  __attribute__((aligned(ALIGN)));
  float* Bx = new float[NN*36+4]  __attribute__((aligned(ALIGN)));
  float* Cx = new float[NN*36+4]  __attribute__((aligned(ALIGN)));
  for (size_t j=0;j<NN*36;++j) Cx[j]=0.;
  //
  // store in element order (i.e. all matrices for a given element are contiguous)
  for (size_t x=0;x<NN;++x) {
    for (size_t i=0;i<36;++i) {
      Ax[i*NN + x] = input[i][x];
      Bx[i*NN + x] = input[i][x]+0.5;
    }
  }

  const clock_t begin = clock();
#pragma vector aligned
  for (size_t x = 0; x < NN; ++x) {
    for (size_t i = 0; i < 6; ++i) {
      for (size_t j = 0; j < 6; ++j) {
        for (size_t k = 0; k < 6; ++k) {
          Cx[ x + (i*6 + j)*NN ] += Ax[ x + (i*6 + k)*NN ] * Bx[ x + (k*6 + j)*NN ];
        }
      }
    }
  }
  const clock_t end = clock();

  // std::cout << "Ax=" << std::endl
  //        << Ax[NN*(0*6+0)] << "\t" << Ax[NN*(0*6+1)] << "\t" << Ax[NN*(0*6+2)] << "\t" << Ax[NN*(0*6+3)] << "\t" << Ax[NN*(0*6+4)] << "\t" << Ax[NN*(0*6+5)] << std::endl
  //        << Ax[NN*(1*6+0)] << "\t" << Ax[NN*(1*6+1)] << "\t" << Ax[NN*(1*6+2)] << "\t" << Ax[NN*(1*6+3)] << "\t" << Ax[NN*(1*6+4)] << "\t" << Ax[NN*(1*6+5)] << std::endl
  //        << Ax[NN*(2*6+0)] << "\t" << Ax[NN*(2*6+1)] << "\t" << Ax[NN*(2*6+2)] << "\t" << Ax[NN*(2*6+3)] << "\t" << Ax[NN*(2*6+4)] << "\t" << Ax[NN*(2*6+5)] << std::endl
  //        << Ax[NN*(3*6+0)] << "\t" << Ax[NN*(3*6+1)] << "\t" << Ax[NN*(3*6+2)] << "\t" << Ax[NN*(3*6+3)] << "\t" << Ax[NN*(3*6+4)] << "\t" << Ax[NN*(3*6+5)] << std::endl
  //        << Ax[NN*(4*6+0)] << "\t" << Ax[NN*(4*6+1)] << "\t" << Ax[NN*(4*6+2)] << "\t" << Ax[NN*(4*6+3)] << "\t" << Ax[NN*(4*6+4)] << "\t" << Ax[NN*(4*6+5)] << std::endl
  //        << Ax[NN*(5*6+0)] << "\t" << Ax[NN*(5*6+1)] << "\t" << Ax[NN*(5*6+2)] << "\t" << Ax[NN*(5*6+3)] << "\t" << Ax[NN*(5*6+4)] << "\t" << Ax[NN*(5*6+5)] << std::endl;
  // std::cout << "Bx=" << std::endl
  //        << Bx[NN*(0*6+0)] << "\t" << Bx[NN*(0*6+1)] << "\t" << Bx[NN*(0*6+2)] << "\t" << Bx[NN*(0*6+3)] << "\t" << Bx[NN*(0*6+4)] << "\t" << Bx[NN*(0*6+5)] << std::endl
  //        << Bx[NN*(1*6+0)] << "\t" << Bx[NN*(1*6+1)] << "\t" << Bx[NN*(1*6+2)] << "\t" << Bx[NN*(1*6+3)] << "\t" << Bx[NN*(1*6+4)] << "\t" << Bx[NN*(1*6+5)] << std::endl
  //        << Bx[NN*(2*6+0)] << "\t" << Bx[NN*(2*6+1)] << "\t" << Bx[NN*(2*6+2)] << "\t" << Bx[NN*(2*6+3)] << "\t" << Bx[NN*(2*6+4)] << "\t" << Bx[NN*(2*6+5)] << std::endl
  //        << Bx[NN*(3*6+0)] << "\t" << Bx[NN*(3*6+1)] << "\t" << Bx[NN*(3*6+2)] << "\t" << Bx[NN*(3*6+3)] << "\t" << Bx[NN*(3*6+4)] << "\t" << Bx[NN*(3*6+5)] << std::endl
  //        << Bx[NN*(4*6+0)] << "\t" << Bx[NN*(4*6+1)] << "\t" << Bx[NN*(4*6+2)] << "\t" << Bx[NN*(4*6+3)] << "\t" << Bx[NN*(4*6+4)] << "\t" << Bx[NN*(4*6+5)] << std::endl
  //        << Bx[NN*(5*6+0)] << "\t" << Bx[NN*(5*6+1)] << "\t" << Bx[NN*(5*6+2)] << "\t" << Bx[NN*(5*6+3)] << "\t" << Bx[NN*(5*6+4)] << "\t" << Bx[NN*(5*6+5)] << std::endl;
  std::cout << "Cx=" << std::endl
            << Cx[NN*(0*6+0)] << "\t" << Cx[NN*(0*6+1)] << "\t" << Cx[NN*(0*6+2)] << "\t" << Cx[NN*(0*6+3)] << "\t" << Cx[NN*(0*6+4)] << "\t" << Cx[NN*(0*6+5)] << std::endl
            << Cx[NN*(1*6+0)] << "\t" << Cx[NN*(1*6+1)] << "\t" << Cx[NN*(1*6+2)] << "\t" << Cx[NN*(1*6+3)] << "\t" << Cx[NN*(1*6+4)] << "\t" << Cx[NN*(1*6+5)] << std::endl
            << Cx[NN*(2*6+0)] << "\t" << Cx[NN*(2*6+1)] << "\t" << Cx[NN*(2*6+2)] << "\t" << Cx[NN*(2*6+3)] << "\t" << Cx[NN*(2*6+4)] << "\t" << Cx[NN*(2*6+5)] << std::endl
            << Cx[NN*(3*6+0)] << "\t" << Cx[NN*(3*6+1)] << "\t" << Cx[NN*(3*6+2)] << "\t" << Cx[NN*(3*6+3)] << "\t" << Cx[NN*(3*6+4)] << "\t" << Cx[NN*(3*6+5)] << std::endl
            << Cx[NN*(4*6+0)] << "\t" << Cx[NN*(4*6+1)] << "\t" << Cx[NN*(4*6+2)] << "\t" << Cx[NN*(4*6+3)] << "\t" << Cx[NN*(4*6+4)] << "\t" << Cx[NN*(4*6+5)] << std::endl
            << Cx[NN*(5*6+0)] << "\t" << Cx[NN*(5*6+1)] << "\t" << Cx[NN*(5*6+2)] << "\t" << Cx[NN*(5*6+3)] << "\t" << Cx[NN*(5*6+4)] << "\t" << Cx[NN*(5*6+5)] << std::endl;
  float time = float(end-begin)/CLOCKS_PER_SEC;
  std::cout << "plainArray_element -- time for NN=" << NN << " multiplications is " << time << " s, i.e. per track [s]=" << time/float(NN) << std::endl;

  delete Ax, Bx, Cx;
}

void test_plainArray_el16mx(float* const* input, const int type) {
  //
  float* Ax = new float[NN*36+4]  __attribute__((aligned(ALIGN)));
  float* Bx = new float[NN*36+4]  __attribute__((aligned(ALIGN)));
  float* Cx = new float[NN*36+4]  __attribute__((aligned(ALIGN)));
  for (size_t j=0;j<NN*36;++j) Cx[j]=0.;
  //
  // store in element order for bunches of 16 matrices (a la matriplex)
  for (size_t x=0;x<NN/16;++x) {
    for (size_t i=0;i<36;++i) {
      for (size_t n=0;n<16;++n) {
        Ax[n + i*16 + 16*36*x] = input[i][n+16*x];
        Bx[n + i*16 + 16*36*x] = input[i][n+16*x]+0.5;
      }
    }
  }

  const clock_t begin = clock();
  if (type==1) {
#pragma vector aligned
    for (size_t x = 0; x < NN/16; ++x) {
      const size_t Nx = x*16*36;
      for (size_t n = 0; n < 16; ++n) {
        Cx[Nx+16* 0+n] = Ax[Nx+16* 0+n]*Bx[Nx+16* 0+n] + Ax[Nx+16* 1+n]*Bx[Nx+16* 6+n] + Ax[Nx+16* 2+n]*Bx[Nx+16*12+n] + Ax[Nx+16* 3+n]*Bx[Nx+16*18+n] + Ax[Nx+16* 4+n]*Bx[Nx+16*24+n] + Ax[Nx+16* 5+n]*Bx[Nx+16*30+n];
        Cx[Nx+16* 1+n] = Ax[Nx+16* 0+n]*Bx[Nx+16* 1+n] + Ax[Nx+16* 1+n]*Bx[Nx+16* 7+n] + Ax[Nx+16* 2+n]*Bx[Nx+16*13+n] + Ax[Nx+16* 3+n]*Bx[Nx+16*19+n] + Ax[Nx+16* 4+n]*Bx[Nx+16*25+n] + Ax[Nx+16* 5+n]*Bx[Nx+16*31+n];
        Cx[Nx+16* 2+n] = Ax[Nx+16* 0+n]*Bx[Nx+16* 2+n] + Ax[Nx+16* 1+n]*Bx[Nx+16* 8+n] + Ax[Nx+16* 2+n]*Bx[Nx+16*14+n] + Ax[Nx+16* 3+n]*Bx[Nx+16*20+n] + Ax[Nx+16* 4+n]*Bx[Nx+16*26+n] + Ax[Nx+16* 5+n]*Bx[Nx+16*32+n];
        Cx[Nx+16* 3+n] = Ax[Nx+16* 0+n]*Bx[Nx+16* 3+n] + Ax[Nx+16* 1+n]*Bx[Nx+16* 9+n] + Ax[Nx+16* 2+n]*Bx[Nx+16*15+n] + Ax[Nx+16* 3+n]*Bx[Nx+16*21+n] + Ax[Nx+16* 4+n]*Bx[Nx+16*27+n] + Ax[Nx+16* 5+n]*Bx[Nx+16*33+n];
        Cx[Nx+16* 4+n] = Ax[Nx+16* 0+n]*Bx[Nx+16* 4+n] + Ax[Nx+16* 1+n]*Bx[Nx+16*10+n] + Ax[Nx+16* 2+n]*Bx[Nx+16*16+n] + Ax[Nx+16* 3+n]*Bx[Nx+16*22+n] + Ax[Nx+16* 4+n]*Bx[Nx+16*28+n] + Ax[Nx+16* 5+n]*Bx[Nx+16*34+n];
        Cx[Nx+16* 5+n] = Ax[Nx+16* 0+n]*Bx[Nx+16* 5+n] + Ax[Nx+16* 1+n]*Bx[Nx+16*11+n] + Ax[Nx+16* 2+n]*Bx[Nx+16*17+n] + Ax[Nx+16* 3+n]*Bx[Nx+16*23+n] + Ax[Nx+16* 4+n]*Bx[Nx+16*29+n] + Ax[Nx+16* 5+n]*Bx[Nx+16*35+n];
        Cx[Nx+16* 6+n] = Ax[Nx+16* 6+n]*Bx[Nx+16* 0+n] + Ax[Nx+16* 7+n]*Bx[Nx+16* 6+n] + Ax[Nx+16* 8+n]*Bx[Nx+16*12+n] + Ax[Nx+16* 9+n]*Bx[Nx+16*18+n] + Ax[Nx+16*10+n]*Bx[Nx+16*24+n] + Ax[Nx+16*11+n]*Bx[Nx+16*30+n];
        Cx[Nx+16* 7+n] = Ax[Nx+16* 6+n]*Bx[Nx+16* 1+n] + Ax[Nx+16* 7+n]*Bx[Nx+16* 7+n] + Ax[Nx+16* 8+n]*Bx[Nx+16*13+n] + Ax[Nx+16* 9+n]*Bx[Nx+16*19+n] + Ax[Nx+16*10+n]*Bx[Nx+16*25+n] + Ax[Nx+16*11+n]*Bx[Nx+16*31+n];
        Cx[Nx+16* 8+n] = Ax[Nx+16* 6+n]*Bx[Nx+16* 2+n] + Ax[Nx+16* 7+n]*Bx[Nx+16* 8+n] + Ax[Nx+16* 8+n]*Bx[Nx+16*14+n] + Ax[Nx+16* 9+n]*Bx[Nx+16*20+n] + Ax[Nx+16*10+n]*Bx[Nx+16*26+n] + Ax[Nx+16*11+n]*Bx[Nx+16*32+n];
        Cx[Nx+16* 9+n] = Ax[Nx+16* 6+n]*Bx[Nx+16* 3+n] + Ax[Nx+16* 7+n]*Bx[Nx+16* 9+n] + Ax[Nx+16* 8+n]*Bx[Nx+16*15+n] + Ax[Nx+16* 9+n]*Bx[Nx+16*21+n] + Ax[Nx+16*10+n]*Bx[Nx+16*27+n] + Ax[Nx+16*11+n]*Bx[Nx+16*33+n];
        Cx[Nx+16*10+n] = Ax[Nx+16* 6+n]*Bx[Nx+16* 4+n] + Ax[Nx+16* 7+n]*Bx[Nx+16*10+n] + Ax[Nx+16* 8+n]*Bx[Nx+16*16+n] + Ax[Nx+16* 9+n]*Bx[Nx+16*22+n] + Ax[Nx+16*10+n]*Bx[Nx+16*28+n] + Ax[Nx+16*11+n]*Bx[Nx+16*34+n];
        Cx[Nx+16*11+n] = Ax[Nx+16* 6+n]*Bx[Nx+16* 5+n] + Ax[Nx+16* 7+n]*Bx[Nx+16*11+n] + Ax[Nx+16* 8+n]*Bx[Nx+16*17+n] + Ax[Nx+16* 9+n]*Bx[Nx+16*23+n] + Ax[Nx+16*10+n]*Bx[Nx+16*29+n] + Ax[Nx+16*11+n]*Bx[Nx+16*35+n];
        Cx[Nx+16*12+n] = Ax[Nx+16*12+n]*Bx[Nx+16* 0+n] + Ax[Nx+16*13+n]*Bx[Nx+16* 6+n] + Ax[Nx+16*14+n]*Bx[Nx+16*12+n] + Ax[Nx+16*15+n]*Bx[Nx+16*18+n] + Ax[Nx+16*16+n]*Bx[Nx+16*24+n] + Ax[Nx+16*17+n]*Bx[Nx+16*30+n];
        Cx[Nx+16*13+n] = Ax[Nx+16*12+n]*Bx[Nx+16* 1+n] + Ax[Nx+16*13+n]*Bx[Nx+16* 7+n] + Ax[Nx+16*14+n]*Bx[Nx+16*13+n] + Ax[Nx+16*15+n]*Bx[Nx+16*19+n] + Ax[Nx+16*16+n]*Bx[Nx+16*25+n] + Ax[Nx+16*17+n]*Bx[Nx+16*31+n];
        Cx[Nx+16*14+n] = Ax[Nx+16*12+n]*Bx[Nx+16* 2+n] + Ax[Nx+16*13+n]*Bx[Nx+16* 8+n] + Ax[Nx+16*14+n]*Bx[Nx+16*14+n] + Ax[Nx+16*15+n]*Bx[Nx+16*20+n] + Ax[Nx+16*16+n]*Bx[Nx+16*26+n] + Ax[Nx+16*17+n]*Bx[Nx+16*32+n];
        Cx[Nx+16*15+n] = Ax[Nx+16*12+n]*Bx[Nx+16* 3+n] + Ax[Nx+16*13+n]*Bx[Nx+16* 9+n] + Ax[Nx+16*14+n]*Bx[Nx+16*15+n] + Ax[Nx+16*15+n]*Bx[Nx+16*21+n] + Ax[Nx+16*16+n]*Bx[Nx+16*27+n] + Ax[Nx+16*17+n]*Bx[Nx+16*33+n];
        Cx[Nx+16*16+n] = Ax[Nx+16*12+n]*Bx[Nx+16* 4+n] + Ax[Nx+16*13+n]*Bx[Nx+16*10+n] + Ax[Nx+16*14+n]*Bx[Nx+16*16+n] + Ax[Nx+16*15+n]*Bx[Nx+16*22+n] + Ax[Nx+16*16+n]*Bx[Nx+16*28+n] + Ax[Nx+16*17+n]*Bx[Nx+16*34+n];
        Cx[Nx+16*17+n] = Ax[Nx+16*12+n]*Bx[Nx+16* 5+n] + Ax[Nx+16*13+n]*Bx[Nx+16*11+n] + Ax[Nx+16*14+n]*Bx[Nx+16*17+n] + Ax[Nx+16*15+n]*Bx[Nx+16*23+n] + Ax[Nx+16*16+n]*Bx[Nx+16*29+n] + Ax[Nx+16*17+n]*Bx[Nx+16*35+n];
        Cx[Nx+16*18+n] = Ax[Nx+16*18+n]*Bx[Nx+16* 0+n] + Ax[Nx+16*19+n]*Bx[Nx+16* 6+n] + Ax[Nx+16*20+n]*Bx[Nx+16*12+n] + Ax[Nx+16*21+n]*Bx[Nx+16*18+n] + Ax[Nx+16*22+n]*Bx[Nx+16*24+n] + Ax[Nx+16*23+n]*Bx[Nx+16*30+n];
        Cx[Nx+16*19+n] = Ax[Nx+16*18+n]*Bx[Nx+16* 1+n] + Ax[Nx+16*19+n]*Bx[Nx+16* 7+n] + Ax[Nx+16*20+n]*Bx[Nx+16*13+n] + Ax[Nx+16*21+n]*Bx[Nx+16*19+n] + Ax[Nx+16*22+n]*Bx[Nx+16*25+n] + Ax[Nx+16*23+n]*Bx[Nx+16*31+n];
        Cx[Nx+16*20+n] = Ax[Nx+16*18+n]*Bx[Nx+16* 2+n] + Ax[Nx+16*19+n]*Bx[Nx+16* 8+n] + Ax[Nx+16*20+n]*Bx[Nx+16*14+n] + Ax[Nx+16*21+n]*Bx[Nx+16*20+n] + Ax[Nx+16*22+n]*Bx[Nx+16*26+n] + Ax[Nx+16*23+n]*Bx[Nx+16*32+n];
        Cx[Nx+16*21+n] = Ax[Nx+16*18+n]*Bx[Nx+16* 3+n] + Ax[Nx+16*19+n]*Bx[Nx+16* 9+n] + Ax[Nx+16*20+n]*Bx[Nx+16*15+n] + Ax[Nx+16*21+n]*Bx[Nx+16*21+n] + Ax[Nx+16*22+n]*Bx[Nx+16*27+n] + Ax[Nx+16*23+n]*Bx[Nx+16*33+n];
        Cx[Nx+16*22+n] = Ax[Nx+16*18+n]*Bx[Nx+16* 4+n] + Ax[Nx+16*19+n]*Bx[Nx+16*10+n] + Ax[Nx+16*20+n]*Bx[Nx+16*16+n] + Ax[Nx+16*21+n]*Bx[Nx+16*22+n] + Ax[Nx+16*22+n]*Bx[Nx+16*28+n] + Ax[Nx+16*23+n]*Bx[Nx+16*34+n];
        Cx[Nx+16*23+n] = Ax[Nx+16*18+n]*Bx[Nx+16* 5+n] + Ax[Nx+16*19+n]*Bx[Nx+16*11+n] + Ax[Nx+16*20+n]*Bx[Nx+16*17+n] + Ax[Nx+16*21+n]*Bx[Nx+16*23+n] + Ax[Nx+16*22+n]*Bx[Nx+16*29+n] + Ax[Nx+16*23+n]*Bx[Nx+16*35+n];
        Cx[Nx+16*24+n] = Ax[Nx+16*24+n]*Bx[Nx+16* 0+n] + Ax[Nx+16*25+n]*Bx[Nx+16* 6+n] + Ax[Nx+16*26+n]*Bx[Nx+16*12+n] + Ax[Nx+16*27+n]*Bx[Nx+16*18+n] + Ax[Nx+16*28+n]*Bx[Nx+16*24+n] + Ax[Nx+16*29+n]*Bx[Nx+16*30+n];
        Cx[Nx+16*25+n] = Ax[Nx+16*24+n]*Bx[Nx+16* 1+n] + Ax[Nx+16*25+n]*Bx[Nx+16* 7+n] + Ax[Nx+16*26+n]*Bx[Nx+16*13+n] + Ax[Nx+16*27+n]*Bx[Nx+16*19+n] + Ax[Nx+16*28+n]*Bx[Nx+16*25+n] + Ax[Nx+16*29+n]*Bx[Nx+16*31+n];
        Cx[Nx+16*26+n] = Ax[Nx+16*24+n]*Bx[Nx+16* 2+n] + Ax[Nx+16*25+n]*Bx[Nx+16* 8+n] + Ax[Nx+16*26+n]*Bx[Nx+16*14+n] + Ax[Nx+16*27+n]*Bx[Nx+16*20+n] + Ax[Nx+16*28+n]*Bx[Nx+16*26+n] + Ax[Nx+16*29+n]*Bx[Nx+16*32+n];
        Cx[Nx+16*27+n] = Ax[Nx+16*24+n]*Bx[Nx+16* 3+n] + Ax[Nx+16*25+n]*Bx[Nx+16* 9+n] + Ax[Nx+16*26+n]*Bx[Nx+16*15+n] + Ax[Nx+16*27+n]*Bx[Nx+16*21+n] + Ax[Nx+16*28+n]*Bx[Nx+16*27+n] + Ax[Nx+16*29+n]*Bx[Nx+16*33+n];
        Cx[Nx+16*28+n] = Ax[Nx+16*24+n]*Bx[Nx+16* 4+n] + Ax[Nx+16*25+n]*Bx[Nx+16*10+n] + Ax[Nx+16*26+n]*Bx[Nx+16*16+n] + Ax[Nx+16*27+n]*Bx[Nx+16*22+n] + Ax[Nx+16*28+n]*Bx[Nx+16*28+n] + Ax[Nx+16*29+n]*Bx[Nx+16*34+n];
        Cx[Nx+16*29+n] = Ax[Nx+16*24+n]*Bx[Nx+16* 5+n] + Ax[Nx+16*25+n]*Bx[Nx+16*11+n] + Ax[Nx+16*26+n]*Bx[Nx+16*17+n] + Ax[Nx+16*27+n]*Bx[Nx+16*23+n] + Ax[Nx+16*28+n]*Bx[Nx+16*29+n] + Ax[Nx+16*29+n]*Bx[Nx+16*35+n];
        Cx[Nx+16*30+n] = Ax[Nx+16*30+n]*Bx[Nx+16* 0+n] + Ax[Nx+16*31+n]*Bx[Nx+16* 6+n] + Ax[Nx+16*32+n]*Bx[Nx+16*12+n] + Ax[Nx+16*33+n]*Bx[Nx+16*18+n] + Ax[Nx+16*34+n]*Bx[Nx+16*24+n] + Ax[Nx+16*35+n]*Bx[Nx+16*30+n];
        Cx[Nx+16*31+n] = Ax[Nx+16*30+n]*Bx[Nx+16* 1+n] + Ax[Nx+16*31+n]*Bx[Nx+16* 7+n] + Ax[Nx+16*32+n]*Bx[Nx+16*13+n] + Ax[Nx+16*33+n]*Bx[Nx+16*19+n] + Ax[Nx+16*34+n]*Bx[Nx+16*25+n] + Ax[Nx+16*35+n]*Bx[Nx+16*31+n];
        Cx[Nx+16*32+n] = Ax[Nx+16*30+n]*Bx[Nx+16* 2+n] + Ax[Nx+16*31+n]*Bx[Nx+16* 8+n] + Ax[Nx+16*32+n]*Bx[Nx+16*14+n] + Ax[Nx+16*33+n]*Bx[Nx+16*20+n] + Ax[Nx+16*34+n]*Bx[Nx+16*26+n] + Ax[Nx+16*35+n]*Bx[Nx+16*32+n];
        Cx[Nx+16*33+n] = Ax[Nx+16*30+n]*Bx[Nx+16* 3+n] + Ax[Nx+16*31+n]*Bx[Nx+16* 9+n] + Ax[Nx+16*32+n]*Bx[Nx+16*15+n] + Ax[Nx+16*33+n]*Bx[Nx+16*21+n] + Ax[Nx+16*34+n]*Bx[Nx+16*27+n] + Ax[Nx+16*35+n]*Bx[Nx+16*33+n];
        Cx[Nx+16*34+n] = Ax[Nx+16*30+n]*Bx[Nx+16* 4+n] + Ax[Nx+16*31+n]*Bx[Nx+16*10+n] + Ax[Nx+16*32+n]*Bx[Nx+16*16+n] + Ax[Nx+16*33+n]*Bx[Nx+16*22+n] + Ax[Nx+16*34+n]*Bx[Nx+16*28+n] + Ax[Nx+16*35+n]*Bx[Nx+16*34+n];
        Cx[Nx+16*35+n] = Ax[Nx+16*30+n]*Bx[Nx+16* 5+n] + Ax[Nx+16*31+n]*Bx[Nx+16*11+n] + Ax[Nx+16*32+n]*Bx[Nx+16*17+n] + Ax[Nx+16*33+n]*Bx[Nx+16*23+n] + Ax[Nx+16*34+n]*Bx[Nx+16*29+n] + Ax[Nx+16*35+n]*Bx[Nx+16*35+n];
      }
    }
  } else {
    for (size_t x = 0; x < NN/16; ++x) {
      const size_t Nx = x*16*36;
      for (size_t i = 0; i < 6; ++i) {
        for (size_t j = 0; j < 6; ++j) {
          for (size_t k = 0; k < 6; ++k) {
            for (size_t n = 0; n < 16; ++n) {
              Cx[ Nx + (i*6 + j)*16 + n ] += Ax[ Nx + (i*6 + k)*16 + n ] * Bx[ Nx + (k*6 + j)*16 + n];
            }
          }
        }
      }
    }
  }
  const clock_t end = clock();

  // std::cout << "Ax=" << std::endl
  //        << Ax[16*(0*6+0)] << "\t" << Ax[16*(0*6+1)] << "\t" << Ax[16*(0*6+2)] << "\t" << Ax[16*(0*6+3)] << "\t" << Ax[16*(0*6+4)] << "\t" << Ax[16*(0*6+5)] << std::endl
  //        << Ax[16*(1*6+0)] << "\t" << Ax[16*(1*6+1)] << "\t" << Ax[16*(1*6+2)] << "\t" << Ax[16*(1*6+3)] << "\t" << Ax[16*(1*6+4)] << "\t" << Ax[16*(1*6+5)] << std::endl
  //        << Ax[16*(2*6+0)] << "\t" << Ax[16*(2*6+1)] << "\t" << Ax[16*(2*6+2)] << "\t" << Ax[16*(2*6+3)] << "\t" << Ax[16*(2*6+4)] << "\t" << Ax[16*(2*6+5)] << std::endl
  //        << Ax[16*(3*6+0)] << "\t" << Ax[16*(3*6+1)] << "\t" << Ax[16*(3*6+2)] << "\t" << Ax[16*(3*6+3)] << "\t" << Ax[16*(3*6+4)] << "\t" << Ax[16*(3*6+5)] << std::endl
  //        << Ax[16*(4*6+0)] << "\t" << Ax[16*(4*6+1)] << "\t" << Ax[16*(4*6+2)] << "\t" << Ax[16*(4*6+3)] << "\t" << Ax[16*(4*6+4)] << "\t" << Ax[16*(4*6+5)] << std::endl
  //        << Ax[16*(5*6+0)] << "\t" << Ax[16*(5*6+1)] << "\t" << Ax[16*(5*6+2)] << "\t" << Ax[16*(5*6+3)] << "\t" << Ax[16*(5*6+4)] << "\t" << Ax[16*(5*6+5)] << std::endl;
  // std::cout << "Bx=" << std::endl
  //        << Bx[16*(0*6+0)] << "\t" << Bx[16*(0*6+1)] << "\t" << Bx[16*(0*6+2)] << "\t" << Bx[16*(0*6+3)] << "\t" << Bx[16*(0*6+4)] << "\t" << Bx[16*(0*6+5)] << std::endl
  //        << Bx[16*(1*6+0)] << "\t" << Bx[16*(1*6+1)] << "\t" << Bx[16*(1*6+2)] << "\t" << Bx[16*(1*6+3)] << "\t" << Bx[16*(1*6+4)] << "\t" << Bx[16*(1*6+5)] << std::endl
  //        << Bx[16*(2*6+0)] << "\t" << Bx[16*(2*6+1)] << "\t" << Bx[16*(2*6+2)] << "\t" << Bx[16*(2*6+3)] << "\t" << Bx[16*(2*6+4)] << "\t" << Bx[16*(2*6+5)] << std::endl
  //        << Bx[16*(3*6+0)] << "\t" << Bx[16*(3*6+1)] << "\t" << Bx[16*(3*6+2)] << "\t" << Bx[16*(3*6+3)] << "\t" << Bx[16*(3*6+4)] << "\t" << Bx[16*(3*6+5)] << std::endl
  //        << Bx[16*(4*6+0)] << "\t" << Bx[16*(4*6+1)] << "\t" << Bx[16*(4*6+2)] << "\t" << Bx[16*(4*6+3)] << "\t" << Bx[16*(4*6+4)] << "\t" << Bx[16*(4*6+5)] << std::endl
  //        << Bx[16*(5*6+0)] << "\t" << Bx[16*(5*6+1)] << "\t" << Bx[16*(5*6+2)] << "\t" << Bx[16*(5*6+3)] << "\t" << Bx[16*(5*6+4)] << "\t" << Bx[16*(5*6+5)] << std::endl;
  std::cout << "Cx=" << std::endl
            << Cx[16*(0*6+0)] << "\t" << Cx[16*(0*6+1)] << "\t" << Cx[16*(0*6+2)] << "\t" << Cx[16*(0*6+3)] << "\t" << Cx[16*(0*6+4)] << "\t" << Cx[16*(0*6+5)] << std::endl
            << Cx[16*(1*6+0)] << "\t" << Cx[16*(1*6+1)] << "\t" << Cx[16*(1*6+2)] << "\t" << Cx[16*(1*6+3)] << "\t" << Cx[16*(1*6+4)] << "\t" << Cx[16*(1*6+5)] << std::endl
            << Cx[16*(2*6+0)] << "\t" << Cx[16*(2*6+1)] << "\t" << Cx[16*(2*6+2)] << "\t" << Cx[16*(2*6+3)] << "\t" << Cx[16*(2*6+4)] << "\t" << Cx[16*(2*6+5)] << std::endl
            << Cx[16*(3*6+0)] << "\t" << Cx[16*(3*6+1)] << "\t" << Cx[16*(3*6+2)] << "\t" << Cx[16*(3*6+3)] << "\t" << Cx[16*(3*6+4)] << "\t" << Cx[16*(3*6+5)] << std::endl
            << Cx[16*(4*6+0)] << "\t" << Cx[16*(4*6+1)] << "\t" << Cx[16*(4*6+2)] << "\t" << Cx[16*(4*6+3)] << "\t" << Cx[16*(4*6+4)] << "\t" << Cx[16*(4*6+5)] << std::endl
            << Cx[16*(5*6+0)] << "\t" << Cx[16*(5*6+1)] << "\t" << Cx[16*(5*6+2)] << "\t" << Cx[16*(5*6+3)] << "\t" << Cx[16*(5*6+4)] << "\t" << Cx[16*(5*6+5)] << std::endl;
  float time = float(end-begin)/CLOCKS_PER_SEC;
  if (type==1) std::cout << "plainArray_el16mx (mplex loop) -- time for NN=" << NN << " multiplications is " << time << " s, i.e. per track [s]=" << time/float(NN) << std::endl;
  else std::cout << "plainArray_el16mx (plain loop) -- time for NN=" << NN << " multiplications is " << time << " s, i.e. per track [s]=" << time/float(NN) << std::endl;
  delete Ax, Bx, Cx;
}
int main(int argc, char* argv[])
{
  float **input = new float*[36];
  for (size_t i = 0; i < 36; i++)  {
    input[i] = new float[NN];
    for (size_t x = 0; x < NN; x++)  
      input[i][x] = (rand()%100)+1;
  }

  std::cout << "done preparing input" << std::endl;

  test_plainArray_matrix(input,0);
  std::cout << std::endl;
  test_plainArray_element(input,0);
  std::cout << std::endl;
  test_plainArray_el16mx(input,0);
  std::cout << std::endl;
  test_plainArray_el16mx(input,1);

  return 0;
}

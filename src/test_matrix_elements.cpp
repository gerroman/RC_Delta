#include "iMe2.h"
#include "iMp2.h"
#include "iMd2.h"
#include "iMpiMe.h"
#include "iMdiMe.h"
#include "iMdiMeSoft.h"
#include "iMdiMeHard.h"
#include "iMdiMp.h"
#include "iMborn2.h"

#include <cstdio>

int main() {

  double nu  =  2.0081047625217883;
  double qe2 = -0.2456103286849673;
  double qp2 = -0.5145880750511680;
  double kK  =  0.2548044526819022;
  double kP  =  1.1100797669282418;
  double kq  =  0.1344888731831004;

  printf(" === check === \n");
  printf("iMe2  check   :%.16f\n", iMe2(nu, qe2, qp2, kK, kP, kq));
  printf("lterm check   :0.1450296502795442\n");
  printf(" === check === \n");
  printf("iMp2  check   :%.16f\n", iMp2(nu, qe2, qp2, kK, kP, kq));
  printf("pterm check   :0.1071087462629029\n");
  printf(" === check === \n");
  printf("iMpiMe  check :%.16f\n", iMpiMe(nu, qe2, qp2, kK, kP, kq));
  printf("lpterm check  :0.1580205013218304\n");
  printf(" === check === \n");
  printf("iMd211 check  :%.16f\n", iMd211(nu, qe2, qp2, kK, kP, kq));
  printf("dterm 11      :0.1763960895328316\n");
  printf(" === check === \n");
  printf("iMd222 check  :%.16f\n", iMd222(nu, qe2, qp2, kK, kP, kq));
  printf("dterm 22      :0.003438344206850841\n");
  printf(" === check === \n");
  printf("iMd212 check  :%.16f\n", iMd212(nu, qe2, qp2, kK, kP, kq));
  printf("dterm 12      :-0.01282306748044593\n");
  printf(" === check === \n");
  printf("iMdiMe1 check :%.16f\n", iMdiMe1(nu, qe2, qp2, kK, kP, kq));
  printf("ldterm1 check :0.0028294752091966\n");
  printf(" === check === \n");
  printf("iMdiMe2 check :%.16f\n", iMdiMe2(nu, qe2, qp2, kK, kP, kq));
  printf("ldterm2 check :-0.0241636369829756\n");
  printf(" === check === \n");
  printf("iMdiMeSoft1   :%.16f\n", iMdiMeSoft1(nu, qe2, qp2, kK, kP, kq));
  printf("ldterms1 check:0.01666898767979055\n");
  printf(" === check === \n");
  printf("iMdiMeSoft2   :%.16f\n", iMdiMeSoft2(nu, qe2, qp2, kK, kP, kq));
  printf("ldterms2 check:-0.00859618431217803\n");
  printf(" === check === \n");
  printf("iMdiMeHard1   :%.16f\n", iMdiMeHard1(nu, qe2, qp2, kK, kP, kq));
  printf("ldtermh1 check:-0.01383951247059386\n");
  printf(" === check === \n");
  printf("iMdiMeHard2   :%.16f\n", iMdiMeHard2(nu, qe2, qp2, kK, kP, kq));
  printf("ldtermh2 check:-0.0155674526707976\n");
  printf(" === check === \n");
  printf("iMdiMp1 check :%.16f\n", iMdiMp1(nu, qe2, qp2, kK, kP, kq));
  printf("iMdiMp2 check :%.16f\n", iMdiMp2(nu, qe2, qp2, kK, kP, kq));
  printf(" === check === \n");
  printf("iMborn2       :%.8f\n", iMborn2(4.47242, -1.51));
  printf("iMborn2 check :0.00282769\n");
  return 0;
}

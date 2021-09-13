CTGSA <- function(x, y){

  soft <- function( cor, alpha ){
    cornew <- max( abs(cor) - alpha, 0 )
    cornew[ cornew < 0 ] <- 0
    sign(cor) * cornew
  }

  mcp <- function( cor, alpha , gamma=2 ){
    if (alpha < cor & cor<= alpha*gamma) {
      mcp <- (gamma/(gamma-1))*(1-alpha/cor)*cor
    } else if (cor > alpha*gamma) {
      mcp <- cor
    } else #(0 < cor & cor <= alpha)
    {mcp <- 0
    }
    return(mcp)
  }

  scad <- function(cor, alpha , gamma=3.7){
    if (0 < cor & cor <= 2*alpha) {
      temp <- max(1-alpha/cor, 0)
      if(temp>=0){
        scad <- temp*cor
      }else{
        scad <- 0
      }
    }else if (2*alpha < cor & cor <= alpha*gamma){
      scad <- (gamma/(gamma-2))*(1-1/gamma-alpha/cor)*cor
    }else{
      scad <- cor
    }
    return(scad)
  }

  p = ncol(x)

  #case,control
  Cx1 = x[which(y==0),]
  Tx1 = x[which(y==1),]

  ####### CTGSA #######
  #alpha(thresholding parameter)
  Cs1 = abs(cor(Cx1));
  Ts1 = abs(cor(Tx1))

  Csv1 = svd(abs(Cs1))$d;
  Tsv1 = svd(abs(Ts1))$d

  #non-coexpression (Ch,Th)
  Ch1 = 1-Csv1[1]/sum(Csv1);
  Th1 = 1-Tsv1[1]/sum(Tsv1)

  H_Cs1 = Cs1
  H_Ts1 = Ts1

  Ctri1 = Cs1[upper.tri(Cs1)];
  Ttri1 = Ts1[upper.tri(Ts1)]

  #non-coexpression(hresholding parameter)
  Ccri1 = sort(Ctri1)[length(Ctri1)*Ch1] ;
  Tcri1 = sort(Ttri1)[length(Ttri1)*Th1]

  #alpha(Ccri,Tcri)thresholding
  #hard
  H_Cs1[Cs1<Ccri1] = 0
  H_Ts1[Ts1<Tcri1] = 0

  Csv1 = svd(abs(H_Cs1))$d;
  Tsv1 = svd(abs(H_Ts1))$d

  #soft
  c_soft1 = soft(H_Cs1,Ccri1);
  t_soft1 = soft(H_Ts1,Tcri1)

  soft_Csv1 = svd(abs(c_soft1))$d;
  soft_Tsv1 = svd(abs(t_soft1))$d

  #mcp
  c_mcp1 = matrix(sapply( H_Cs1, function(z) mcp(z, Ccri1, gamma=2.0) ), p, p, byrow=FALSE )
  t_mcp1 = matrix(sapply( H_Ts1, function(z) mcp(z, Tcri1, gamma=2.0) ), p, p, byrow=FALSE )

  mcp_Csv1 = svd(abs(c_mcp1))$d;
  mcp_Tsv1 = svd(abs(t_mcp1))$d

  #scad
  c_scad1 = matrix(sapply( H_Cs1, function(z) scad(z, Ccri1, gamma=2.0) ), p, p, byrow=FALSE )
  t_scad1 = matrix(sapply( H_Ts1, function(z) scad(z, Tcri1, gamma=2.0) ), p, p, byrow=FALSE )

  scad_Csv1 = svd(abs(c_scad1))$d;
  scad_Tsv1 = svd(abs(t_scad1))$d

  T_soft1=abs(soft_Tsv1[1] - soft_Csv1[1])
  T_mcp1=abs(mcp_Tsv1[1] - mcp_Csv1[1])
  T_scad1=abs(scad_Tsv1[1] - scad_Csv1[1])
  T_TH1 = abs(Tsv1[1] - Csv1[1])

  list(stat=list(hard=T_TH1,
                 soft=T_soft1,
                 mcp=T_mcp1,
                 scad=T_scad1),
       params=list(x=x, y=y))
}

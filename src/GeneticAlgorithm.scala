import scala.math.{cos, sin, abs, Pi}

class System (sizec: Int) {
  val size = sizec
  
  def calculateFit(id: Int) : Float = {
    0.0f
  }
}

class FourierFunctionFit (
    sizec: Int,
    ncoeffsc: Int, 
    x: List[Float], 
    fMatch: Float => Float) extends System(sizec) {
  val ncoeffs = ncoeffsc
  val fMatchValues = List.tabulate(x.size)(n => fMatch(x(n)))
  
  //sines(m)(n) = sines(n*x(m)) = sine weight for coefficient n at m'th x. 
  //Same for cosine.
  val sines = List.tabulate(x.size)(m => List.tabulate(ncoeffs)(n => sin((n+1)*x(m)).toFloat))
  val cosines = List.tabulate(x.size)(m => List.tabulate(ncoeffs)(n => cos((n+1)*x(m)).toFloat))

  //constant shifts
  val a0 = List.fill(size)(0.0f)
  
  //coeffs(i)(n) = n'th coefficient set for individual i 
  val sineCoeffs = List.fill(size)(Array.fill(ncoeffs)(0.0f))
  val cosineCoeffs = List.fill(size)(Array.fill(ncoeffs)(0.0f))

  
  def fourierElement(coeffs: Array[Float], trigioms: List[Float]) : Float = {
      (coeffs, trigioms).zipped.map(_ * _).sum
  }
  
  def calculateFourierFunction(id: Int) : List[Float] = {
     (cosines, sines).zipped.map((c, s) => a0(id) + fourierElement(cosineCoeffs(id), c) + fourierElement(sineCoeffs(id), s))
  }
  
  final override def calculateFit (id: Int) : Float = {
    (fMatchValues, calculateFourierFunction(id)).zipped.map((f, g) => (f-g)*(f-g)).sum 
  }
  
}

class Individual (idc: Int) {
  val id = idc
}

class Population (sizec: Int) {
  val size: Int = sizec
}


object GeneticAlgorithm {
  
  def main(args: Array[String]): Unit = {
    val size = 10
    val ncoeffs = 2
    
    val xn = 100
    val xMax = 2.0f*Pi.toFloat
    
    val x = List.tabulate(xn)(n => xMax*(n/((xn - 1).toFloat)))
    
    //testcase with f(x) = sin(x) + 5cos(x)
    val fff = new FourierFunctionFit(size, ncoeffs, x, (x: Float) => sin(x).toFloat + 5.0f*cos(2*x).toFloat)
    fff.sineCoeffs(0)(0) = 1.0f
    fff.cosineCoeffs(0)(1) = 5.0f
    
    val m = fff.calculateFit(0)
    
    var i = 0
    for (i <- 0 to xn - 1) {
     
      println(fff.fourierElement(fff.sineCoeffs(0), fff.sines(i)), fff.fMatchValues(i))
     
    }
    
    println(m)
    
  }
  
}
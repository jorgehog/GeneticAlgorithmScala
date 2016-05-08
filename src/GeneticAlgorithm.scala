import scala.math.{cos, sin, abs, Pi, floor}
import java.io.{PrintWriter, File}
import scala.util.Random


abstract class Population (sizec: Int) {
  val size = sizec
  
  def calculateFit(id: Int) : Float 
  
  def mutateIndividual(id: Int) : Unit
  
  def generateOffspring(idFirstParent: Int, idSecondParent: Int, 
        idChild: Int) : Unit 
  
  def dumpFile(n: Int, id: Int) : Unit = {}

  
  def inplaceSort(ids: Array[Int], fits: Array[Float]) : Unit = {
    ids.sortWith(fits(_) < fits(_)).zipWithIndex.foreach {case(sortedId, idx) => ids(idx) = sortedId}
  }
  
  def evolve(dumpInterval: Int = 10, nMax: Int = 10000, convErr: Float = 1e-3f) : Int = {
    
    val sortedIds = Array.range(0, size)
    val fits = Array.tabulate(size)(n => calculateFit(n))
    var n = 0
    
    inplaceSort(sortedIds, fits)
    
    while (fits(sortedIds(0)) > convErr & n < nMax) {
    
      //we let the winner generate offsprings with the rest of the
      //half winning population. The last half is replaced by these.
      for (i <- Range(1, size/2)) {
        val child = sortedIds(size - i)
        generateOffspring(sortedIds(0), sortedIds(i), child)
        mutateIndividual(child)
        fits(child) = calculateFit(child)
      }
      
      inplaceSort(sortedIds, fits)
      
      if (n % dumpInterval == 0) {
        println("Generation " + n.toString)
        dumpFile(n/dumpInterval, sortedIds(0))
      }
      
      n += 1
      
    } 
      
    println("n generations: " + n.toString() + " Winner fit: " + fits(sortedIds(0)).toString())
    
    return sortedIds(0)
    
  }
  
}

class FourierSeriesFit (
    sizec: Int,
    ncoeffsc: Int, 
    x: List[Float], 
    targetValuesc: List[Float]) extends Population(sizec) {
  val ncoeffs = ncoeffsc
  val targetValues = targetValuesc
  
  //sines(m)(n) = sines(n*x(m)) = sine weight for coefficient n at m'th x. 
  //Same for cosine.
  val sines = List.tabulate(x.size)(m => List.tabulate(ncoeffs)(n => sin((n+1)*x(m)).toFloat))
  val cosines = List.tabulate(x.size)(m => List.tabulate(ncoeffs)(n => cos((n+1)*x(m)).toFloat))

  //constant shifts
  val a0 = Array.tabulate(size)(n => Random.nextFloat())
  
  //coeffs(i)(n) = n'th coefficient set for individual i 
  val sineCoeffs = List.fill(size)(Array.tabulate(ncoeffs)(n => Random.nextFloat()))
  val cosineCoeffs = List.fill(size)(Array.tabulate(ncoeffs)(n => Random.nextFloat()))

  
  /*
   * Member functions:
   */
  def fourierElement(coeffs: Array[Float], trigioms: List[Float]) : Float = {
    (coeffs, trigioms).zipped.map(_ * _).sum
  }
  
  def calculateFourierFunction(id: Int) : List[Float] = {
    (cosines, sines).zipped.map((c, s) => a0(id) + fourierElement(cosineCoeffs(id), c) + fourierElement(sineCoeffs(id), s))
  }
  
  def mutateSingle(coeff: Float, changeFactor: Float) : Float = {
    changeFactor*coeff
  }
  
  //50-50 mix of each parent
  def mixParents(coeffs: List[Array[Float]], idFirstParent: Int, 
                   idSecondParent: Int, idChild: Int) : Unit = {
    (coeffs(idFirstParent), coeffs(idSecondParent), Range(0, ncoeffs))
      .zipped.foreach((s1, s2, i) => coeffs(idChild)(i) = (s1+s2)/2)
  }
  
  
  /*
   * System Interface:
   */
  final override def calculateFit (id: Int) : Float = {
    (targetValues, calculateFourierFunction(id)).zipped.map((f, g) => (f-g)*(f-g)).sum/targetValues.size
  }
  
  final override def mutateIndividual (id: Int) : Unit = {
      
    //0 => a0, 1-ncoeff => sine coeff, > ncoeff => cosine coeff
    val coefficient = floor((2*ncoeffs+1)*Random.nextFloat()).toInt
    
    val changeFactor = 1.0f + Random.nextGaussian().toFloat
      
    if (coefficient == 0)
    {
      a0(id) = mutateSingle(a0(id), changeFactor)
    }
    else if (coefficient > ncoeffs) {
      val relCoeff = coefficient - ncoeffs - 1
      val coeff = cosineCoeffs(id)(relCoeff)
      cosineCoeffs(id)(relCoeff) = mutateSingle(coeff, changeFactor)
    }
    else {
      val relCoeff = coefficient - 1
      val coeff = sineCoeffs(id)(relCoeff)
      sineCoeffs(id)(relCoeff) = mutateSingle(coeff, changeFactor)
    }
  }
  
  final override def generateOffspring (idParent1: Int, idParent2: Int, idChild: Int) : Unit = {
    a0(idChild) = (a0(idParent1) + a0(idParent2))/2
    
    mixParents(sineCoeffs, idParent1, idParent2, idChild)
    mixParents(cosineCoeffs, idParent1, idParent2, idChild)
  }
  
  final override def dumpFile(n: Int, id: Int) : Unit = {
    val writer = new PrintWriter(new File(f"/tmp/genetic$n%d.dat")) 

    writer.write(a0(id).toString() + "\n")
    (cosineCoeffs(id), sineCoeffs(id)).zipped.foreach{(c, s) => writer.write(f"$c%f $s%f\n")}
    
    writer.close()
  }

}

object GeneticAlgorithm {
  
  def dumpTargetData(x: List[Float], targetData: List[Float]) = {
    val writer = new PrintWriter(new File(f"/tmp/target_data.dat")) 

    (x, targetData).zipped.foreach{(xi, fi) => writer.write(f"$xi%f $fi%f\n")}
    
    writer.close()
  }
  
  def main(args: Array[String]): Unit = {
    val size = 100
    val ncoeffs = 10
   
    val xn = 1000
   
    //testcase with f(x) = {0 if x<0.5 else 1
    def targetFunction = (x: Float) => if (x > 0.0f) 1.0f else 0.0f
    val xMin = -0.5f
    val xMax = 0.5f
    
    //testcase with f(x) = 0.5 + sin(x) + 5cos(x) + 2sin(5x)
    //def targetFunction = (x: Float) => 0.5f + sin(x).toFloat + 5.0f*cos(2*x).toFloat + 2*sin(5*x).toFloat 
    //val xMin = 0.0f
    //val xMax = 2.0f*Pi.toFloat
    
    val x = List.tabulate(xn)(n => xMin + (xMax-xMin)*(n/((xn - 1).toFloat))) 
    val targetValues = List.tabulate(x.size)(n => targetFunction(x(n)))
    
    dumpTargetData(x, targetValues)
    
    val fff = new FourierSeriesFit(size, ncoeffs, x, targetValues)
    
    val winner = fff.evolve(1)
    
    println(fff.a0(winner))
    (fff.sineCoeffs(winner), fff.cosineCoeffs(winner)).zipped.foreach((s, c) => println(s, c))
  }
}

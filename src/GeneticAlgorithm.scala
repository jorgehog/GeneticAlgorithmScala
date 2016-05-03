import scala.math.{cos, sin, abs, Pi, floor}
import java.io.{PrintWriter, File}

abstract class System (sizec: Int) {
  val size = sizec
  
  def calculateFit(id: Int) : Float 
  
  def mutateIndividual(id: Int) : Unit
  
  def generateOffspring(idFirstParent: Int, idSecondParent: Int, idChild: Int) : Unit 
  
  def dumpFile(n: Int, id: Int) : Unit = {}
  
}

class FourierSeriesFit (
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
  val sineCoeffs = List.fill(size)(Array.tabulate(ncoeffs)(n => scala.util.Random.nextFloat()))
  val cosineCoeffs = List.fill(size)(Array.tabulate(ncoeffs)(n => scala.util.Random.nextFloat()))

  
  /*
   * Member functions:
   */
  def fourierElement(coeffs: Array[Float], trigioms: List[Float]) : Float = {
      (coeffs, trigioms).zipped.map(_ * _).sum
  }
  
  def calculateFourierFunction(id: Int) : List[Float] = {
     (cosines, sines).zipped.map((c, s) => a0(id) + fourierElement(cosineCoeffs(id), c) + fourierElement(sineCoeffs(id), s))
  }
  
  def mutateSingle(coeff: Float, phaseShift: Float, amplitudeShift: Float) : Float = {
    //amplitudeShift*coeff
    amplitudeShift*(coeff + phaseShift)
  }
  
  //50-50 mix of each parent
  def generateOffspringCoefficient(coeffs: List[Array[Float]], idFirstParent: Int, 
                                   idSecondParent: Int, idChild: Int) : Unit = {
    (coeffs(idFirstParent), coeffs(idSecondParent), (0 to ncoeffs))
      .zipped.foreach((s1: Float, s2: Float, i: Int) => coeffs(idChild)(i) = (s1+s2)/2)
  }
  
  
  /*
   * System Interface:
   */
  final override def calculateFit (id: Int) : Float = {
    (fMatchValues, calculateFourierFunction(id)).zipped.map((f, g) => (f-g)*(f-g)).sum/fMatchValues.size
  }
  
  final override def mutateIndividual (id: Int) : Unit = {
      
    val coefficient = floor(ncoeffs*scala.util.Random.nextFloat()).toInt
    val phaseShift = scala.util.Random.nextGaussian().toFloat
    val amplitudeShift = 1.0f + scala.util.Random.nextGaussian().toFloat
      
    if (scala.util.Random.nextFloat() < 0.5) {
      sineCoeffs(id)(coefficient) = mutateSingle(sineCoeffs(id)(coefficient), phaseShift, amplitudeShift)
      }
    else {
      cosineCoeffs(id)(coefficient) = mutateSingle(cosineCoeffs(id)(coefficient), phaseShift, amplitudeShift)
    }
  }
  
  final override def generateOffspring (idFirstParent: Int, idSecondParent: Int, idChild: Int) : Unit = {
    generateOffspringCoefficient(sineCoeffs, idFirstParent, idSecondParent, idChild)
    generateOffspringCoefficient(cosineCoeffs, idFirstParent, idSecondParent, idChild)
  }
  
  final override def dumpFile(n: Int, id: Int) : Unit = {
    println("dumping...", n)
    val writer = new PrintWriter(new File(f"/tmp/genetic$n%d.dat")) 

    (cosineCoeffs(id), sineCoeffs(id)).zipped.foreach{(c, s) => writer.write(f"$c%f $s%f\n")}
    
    writer.close()
  }

}

object GeneticAlgorithm {

  def inplaceSort(ids: Array[Int], fits: Array[Float]) : Unit = {
    ids.sortWith(fits(_) < fits(_)).zipWithIndex.foreach {x => ids(x._2) = x._1}
  }
  
  def evolve(system: System, dumpInterval: 
      Int = 10, nMax: Int = 10000, convErr: Float = 1e-3f) : Int = {
    
    val sortedIds = Array.range(0, system.size)
    val fits = Array.tabulate(system.size)(n => system.calculateFit(n))
    var n = 0
    
    inplaceSort(sortedIds, fits)
    
    do {
    
      //we let the winner generate offsprings with the rest of the
      //half winning population. The last half is replaced by these.
      for (i <- Range(1, system.size/2)) {
        val child = sortedIds(system.size - i)
        system.generateOffspring(sortedIds(0), sortedIds(i), child)
        system.mutateIndividual(child)
        fits(child) = system.calculateFit(child)
      }
      
      inplaceSort(sortedIds, fits)
    
      println("trying...", n)
      
      if (n % dumpInterval == 0) {
        system.dumpFile(n/dumpInterval, sortedIds(0))
      }
      
      n += 1
      
    } while (fits(sortedIds(0)) > convErr & n < nMax)
      
    println("n generations: ", n, " Winner fit: ", fits(sortedIds(0)))
    
    return sortedIds(0)
    
  }
  
  def main(args: Array[String]): Unit = {
    val size = 100
    val ncoeffs = 20
    
    val xn = 1000
    //val xMax = 2.0f*Pi.toFloat
    val xMax = 0.5f
    val xMin = -0.5f
    
    val x = List.tabulate(xn)(n => xMin + (xMax-xMin)*(n/((xn - 1).toFloat)))
    
    //testcase with f(x) = sin(x) + 5cos(x)
    //def targetFunction = (x: Float) => sin(x).toFloat + 5.0f*cos(2*x).toFloat + 2*sin(5*x).toFloat
    def targetFunction = (x: Float) => if (x > 0.0f) 1.0f else 0.0f
    val fff = new FourierSeriesFit(size, ncoeffs, x, targetFunction)
    
    val winner = evolve(fff)
    
    (fff.sineCoeffs(winner), fff.cosineCoeffs(winner)).zipped.foreach((s, c) => println(s, c))
  }
}

using System;
using System.IO;
using System.Collections.Generic;

using Pointer = System.Int32;
using static Kakeya.Utility;

namespace Kakeya
{
    public class Program
    {
        public static void Main()
        {
            Random rand = new Random();
            while (true)
            {
                try
                {
                    Console.WriteLine("Prime number:");
                    int p = int.Parse(Console.ReadLine());
                    if (!IsPrime(p))
                    {
                        throw new FormatException("The current given should have been a prime number.");
                    }
                    Console.WriteLine("Power:");
                    int m = int.Parse(Console.ReadLine());
                    if (m <= 0)
                    {
                        throw new FormatException("The current given should have been a positive integer.");
                    }
                    Console.WriteLine("Dimension:");
                    int n = int.Parse(Console.ReadLine());
                    if (n <= 0)
                    {
                        throw new FormatException("The current given should have been a positive integer.");
                    }
                    Console.WriteLine("Number of trials:");
                    int trials = int.Parse(Console.ReadLine());
                    Polynomial f = new Polynomial(m);
                    if (m != 1)
                    {
                        Console.WriteLine("Searching for a polynomial of degree " + m + " irreducible over F_" + p + ".");
                        f = IrreducibleFunction(p, m);
                        Console.WriteLine("Found irreducible polynomial: " + f.Show());
                    }

                    int count = int.MaxValue;
                    int trial = 0;
                    List<Kakeya> Ks = new List<Kakeya>();
                    for (int t = 1; t <= Math.Max(1, trials); t++)
                    {
                        Console.WriteLine("Trial " + t + ":");
                        Kakeya K = new Kakeya(p, m, n, f, rand);
                        int currentCount = K.Print(false);
                        if (currentCount < count)
                        {
                            Ks.Clear();
                            Ks.Add(K);
                            count = currentCount;
                            trial = t;
                        }
                    }
                    Console.WriteLine("\nBest result: (Trial " + trial + ")");
                    foreach (Kakeya K in Ks)
                    {
                        K.Print(true);
                    }

                }
                catch (FormatException e)
                {
                    Console.WriteLine(e.Message);
                    Console.WriteLine("Terminating program...");
                    break;
                }
            }
        }
    }

    public class Kakeya
    {
        int p;
        int q;
        int n;
        Polynomial f;
        Random rand;
        bool[] grid;

        public Kakeya(int p, int m, int n, Polynomial f, Random rand)
        {
            this.p = p;
            this.q = (int)Math.Pow(p, m);
            this.n = n;
            this.f = f;
            this.rand = rand;
            grid = new bool[(int)Math.Pow(q, n)];

            List<int> V = GetList(1, (int)Math.Pow(q, n) - 1); // Think of V as the set of (non-zero) directions in F_q^n to be covered in the Kakeya set, but presented as integers.

            while (V.Count > 0)
            {
                // Take a direction in V:
                int j = GetElement(V, rand);
                Container<Tuple> v = IntToTuples(j, p, q, n);

                // Take a point P such that its corresponding line with direction v has the maximal amount of intersections with the grid.
                // Is slower but generally results in significantly smaller sets than the naive approach.
                Container<Tuple> t = GreedyPoint(v);
                for (int i = 0; i < q; i++)
                {
                    Tuple x = IntToTuple(i, p); // Think of x as a scalar in F_q.

                    Container<Tuple> u = v.Times(x, p).Mod(f, p);
                    V.Remove(u.ToInt(p, q));

                    Container<Tuple> s = t.Plus(u, p);

                    int k = s.ToInt(p, q);

                    grid[k] = true;
                }
            }
        }

        // Returns a point P for which the line with direction v has maximal intersections with the grid.
        private Container<Tuple> GreedyPoint(Container<Tuple> v)
        {
            List<Container<Tuple>> greedyPoints = new List<Container<Tuple>>();
            int greedyCount = 0;

            int total = (int)Math.Pow(q, n);
            bool[] skippable = new bool[total];
            for (int i = 0; i < total; i++)
            {
                if (!skippable[i])
                {
                    Container<Tuple> u = IntToTuples(i, p, q, n);
                    int count = 0;
                    for (int j = 0; j < q; j++)
                    {
                        Tuple t = IntToTuple(j, p);
                        Container<Tuple> c = u.Plus(v.Times(t, p).Mod(f, p), p);
                        int k = c.ToInt(p, q);

                        if (grid[k])
                        {
                            count++;
                        }

                        // As the point corresponding to the integer k lies on the same line as the point corresponding to i, we may skip it in the future.
                        skippable[k] = true;
                    }
                    if (count >= greedyCount)
                    {
                        if (count > greedyCount)
                        {
                            greedyPoints.Clear();
                            greedyCount = count;
                        }
                        greedyPoints.Add(u);
                    }

                }
            }

            return GetElement(greedyPoints, rand);
        }

        public int Print(bool save)
        {
            string res;
            int count = 0;

            Write("q = " + q + ", n = " + n + ":\n", save, false);
            for (int layer = 0; layer < Math.Max(1, (int)Math.Pow(q, n - 2)); layer++)
            {
                int oldCount = count;
                res = "";
                if (n > 2)
                {
                    res += "Layer " + IntToTuple(layer, q, n - 2).Show() + ":\n";
                }
                int total = (int)Math.Pow(q, Math.Min(2, n));
                for (int i = layer * total; i < (layer + 1) * total; i++)
                {
                    if (grid[i])
                    {
                        res += "+";
                        count++;
                    }
                    else
                    {
                        res += " ";
                    }
                    if (i % q == q - 1)
                    {
                        res += "\n";
                    }
                }
                if (oldCount != count)
                {
                    Write(res, save && n < 5, n < 5);
                }
            }
            res = "Count: " + count.ToString() + "\n";
            Write(res, save, true);

            return count;
        }

        private static void Write(string s, bool save, bool print)
        {
            if (print)
            {
                Console.WriteLine(s);
            }
            if (save)
            {
                string binDirectory = Directory.GetParent(Environment.CurrentDirectory).Parent.FullName;
                using (var file = new StreamWriter(binDirectory + "/results.txt", true))
                {
                    file.WriteLine(s);
                    file.Close();
                }
            }
        }
    }

    public class Utility
    {
        // Find integer x such that a * x = b (mod p). We assume p does not divide a.
        public static int FindInverse(int a, int b, int p)
        {
            // We first use the extended Euclidean algorithm to find an integer x satisfying a * x = 1 (mod p).
            // See https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm for more details.
            int prev_r = a, r = p;
            int prev_s = 1, s = 0;

            while (r != 0)
            {
                int q = prev_r / r;
                (prev_r, r) = (r, prev_r - q * r);
                (prev_s, s) = (s, prev_s - q * s);
            }

            // Now we have that a * prev_s = 1 (mod p). Multiplying by b and reducing modulo p, we are done.
            return prev_s * b % p;
        }

        public static Polynomial IrreducibleFunction(int p, int deg)
        {
            Polynomial f = new Polynomial(0);

            int lower = (int)Math.Pow(p, deg); // Corresponds to the first polynomial of degree deg (which is monic).
            int upper = lower << 1; // I.e., upper = 2 * lower (corresponds to the first non-monic polynomial of degree deg).
            for (int n = lower; n < upper; n++)
            {
                f = IntToPolynomial(n, p);
                if (f.Irreducible(p))
                {
                    break;
                }
            }
            return f;
        }

        public static List<int> GetList(int min, int max)
        {
            List<int> L = new List<int>();
            for (int i = min; i <= max; i++)
            {
                L.Add(i);
            }
            return L;
        }

        public static T GetElement<T>(List<T> L, Random rand)
        {
            int index = rand.Next(L.Count);
            T result = L[index];
            return result;
        }

        // From wikipedia: https://en.wikipedia.org/wiki/Primality_test
        public static bool IsPrime(int n)
        {
            if (n == 2 || n == 3)
            {
                return true;
            }
            if (n <= 1 || n % 2 == 0 || n % 3 == 0)
            {
                return false;
            }
            for (int i = 5; i * i <= n; i += 6)
            {
                if (n % i == 0 || n % (i + 2) == 0)
                {
                    return false;
                }
            }
            return true;
        }

        public static Tuple IntToTuple(int n, int b, int size = 1)
        {
            return IntToContainer(n, b, x => new Tuple(x), size);
        }

        public static Polynomial IntToPolynomial(int n, int b)
        {
            return IntToContainer(n, b, x => new Polynomial(x));
        }

        private static C IntToContainer<C>(int n, int b, Func<int, C> f, int size = 1) where C : Container
        {
            // Idea: store the integer n as a container with its entries given by writing n in base b.
            // The number of entries of the container is equal to the smallest integer m such that b^m >= n.
            int k = (int)Math.Pow(b, size);
            while (n >= k)
            {
                k *= b;
                size++;
            }

            C c = f(size);
            for (int i = 0; i < size; i++)
            {
                c.Add(n % b);
                n /= b;
            }
            return c;
        }

        public static Container<Tuple> IntToTuples(int n, int mod1, int mod2, int d)
        {
            Container<Tuple> C = new Container<Tuple>(d);
            for (int i = 0; i < d; i++)
            {
                Tuple t = IntToTuple(n % mod2, mod1);
                C.Add(t);
                n /= mod2;
            }
            return C;
        }
    }

    public class Container<T> where T : Container
    {
        public T[] c;
        private Pointer p;

        public Container(int n, int startIndex = 0)
        {
            c = new T[n];
            p = startIndex;
        }

        public void Add(T n)
        {
            if (p < c.Length)
            {
                c[p++] = n;
            }
            else
            {
                throw new IndexOutOfRangeException("The container is full.");
            }
        }

        public Container<U> Map<U>(Func<T, U> f) where U : Container
        {
            Container<U> C = new Container<U>(p);
            for (int i = 0; i < p; i++)
            {
                C.Add(f(c[i]));
            }
            return C;
        }

        public Container<T> Plus(Container<T> c1, int mod = int.MaxValue)
        {
            if (p < c1.p)
            {
                return c1.Plus(this, mod);
            }

            Container<T> C = new Container<T>(p);
            for (int i = 0; i < c1.p; i++)
            {
                C.Add((T)c[i].Plus(c1.c[i], mod));
            }
            for (int i = c1.p; i < p; i++)
            {
                C.Add(c[i]);
            }
            return C;
        }

        public Container<T> Times(Container c1, int mod = int.MaxValue)
        {
            return Map(x => (T)x.Times(c1, mod));
        }

        public Container<T> Mod(Container c1, int mod = int.MaxValue)
        {
            return Map(x => (T)x.Mod(c1, mod));
        }

        public int ToInt(int mod1, int mod2)
        {
            int sum = 0;
            for (int i = p - 1; i >= 0; i--)
            {
                sum *= mod2;
                sum += c[i].ToInt(mod1);
            }
            return sum;
        }
    }

    public class Container
    {
        public int[] c;
        protected Pointer p; // Represents pointer to first unused index of c.

        public Container(int n, int startIndex = 0)
        {
            c = new int[n];
            p = startIndex;
        }

        public void Add(int n)
        {
            if (p < c.Length)
            {
                c[p++] = n;
            }
            else
            {
                throw new IndexOutOfRangeException("The container is full.");
            }
        }

        public int Last()
        {
            if (p > 0)
            {
                return c[p - 1];
            }
            else
            {
                throw new IndexOutOfRangeException("The container is empty.");
            }
        }

        protected int Size()
        {
            return c.Length;
        }

        public void AdjustCoefficients(int b)
        {
            for (int i = 0; i < p; i++)
            {
                if (c[i] < 0)
                {
                    c[i] += b;
                }
            }
        }

        public Container Plus(Container c1, int mod = int.MaxValue)
        {
            if (p < c1.p)
            {
                return c1.Plus(this, mod);
            }

            Container C = NewContainer(p);
            for (int i = 0; i < c1.p; i++)
            {
                C.Add((c[i] + c1.c[i]) % mod);
            }
            for (int i = c1.p; i < p; i++)
            {
                C.Add(c[i] % mod);
            }
            C.UpdatePointer();
            return C;
        }

        public Container Times(Container c1, int mod = int.MaxValue)
        {
            if (Size() < c1.Size())
            {
                return c1.Times(this, mod);
            }

            Container C = NewContainer(p + c1.p - 1);
            for (int n = 0; n < C.Size(); n++)
            {
                int sum = 0;
                for (int j = Math.Max(0, n - (p - 1)); j <= Math.Min(n, c1.p - 1); j++)
                {
                    sum += c[n - j] * c1.c[j] % mod;
                }
                C.Add(sum % mod);
            }
            C.UpdatePointer();
            return C;
        }

        public Container Mod(Container f, int mod = int.MaxValue)
        {
            Container C = this;
            if (f.Size() > 1)
            {
                while (C.p >= f.Size())
                {
                    Container fCopy = f.Copy(C.p, C.p - f.Size());
                    int lastC = C.Last();
                    int lastF = fCopy.Last();
                    int value = FindInverse(lastF, lastC, mod);
                    fCopy.Times(-value, mod);
                    C = C.Plus(fCopy, mod);
                }
            }
            C.UpdatePointer();
            return C;
        }

        private void UpdatePointer()
        {
            for (int i = p - 1; i >= 0; i--)
            {
                if (c[i] != 0 || i == 0)
                {
                    p = i + 1;
                    break;
                }
            }
        }

        private void Times(int scalar, int mod = int.MaxValue)
        {
            int oldPointer = p;
            for (int i = 0; i < oldPointer; i++)
            {
                c[i] *= scalar;
                c[i] %= mod;
            }
        }

        public int ToInt(int b)
        {
            // First we ensure that the coefficients of the container are between 0 and p - 1 (inclusive).
            AdjustCoefficients(b);

            // Now we convert the container into an integer simply by thinking of the container as a number written in base b.
            int sum = 0;
            for (int i = p - 1; i >= 0; i--)
            {
                sum *= b;
                sum += c[i];
            }
            return sum;
        }

        public virtual string Show()
        {
            string s = "";
            if (p != 1)
            {
                s += "(";
            }
            for (int i = p - 1; i >= 0; i--)
            {
                s += c[i];
                if (i != 0)
                {
                    s += ", ";
                }
            }
            if (p != 1)
            {
                s += ")";
            }
            return s;
        }

        private Container Copy(int length, int startIndex = 0)
        {
            Container C = NewContainer(length, startIndex);
            for (int i = 0; i < p; i++)
            {
                C.Add(c[i]);
            }
            return C;
        }

        public Container NewContainer(int n, int startIndex = 0)
        {
            if (this is Tuple)
            {
                return new Tuple(n, startIndex);
            }
            else if (this is Polynomial)
            {
                return new Polynomial(n, startIndex);
            }
            else return new Container(n, startIndex);
        }
    }

    public class Tuple : Container
    {
        public Tuple(int n, int startIndex = 0) : base(n, startIndex)
        {

        }
    }

    public class Polynomial : Container
    {
        // A polynomial can be saved as a container.
        // For example: 2 + 3X^2 + X^3 <-> {2, 0, 3, 1}
        public Polynomial(int n, int startIndex = 0) : base(n, startIndex)
        {

        }

        public bool IsMonic()
        {
            return Last() == 1;
        }

        public bool IsZero()
        {
            for (int i = 0; i < p; i++)
            {
                if (c[i] != 0)
                {
                    return false;
                }
            }
            return true;
        }

        public int Degree()
        {
            return p - 1;
        }

        public bool Irreducible(int mod)
        {
            // If the degree of f is 2 or 3, then f is reducible precisely when f has a linear factor which is easy to check.
            for (int x = 0; x < mod; x++)
            {
                int sum = 0;
                int currentPower = 1;
                for (int j = 0; j <= Degree(); j++)
                {
                    sum = (sum + c[j] * currentPower) % mod;
                    currentPower = currentPower * x % mod;
                }
                if (sum == 0)
                {
                    return false;
                }
            }
            // If the degree of f is at least 4, then the previous does not guarantee that f is irreducible.
            // In that case we still need to check that f is not the product of two polynomials of degree 2 or higher. We check this by exhaustion,
            // by checking that any polynomial of degree 2 or higher is not a factor of f. Note that we do not have to check above degree higher
            // than deg(f)/2.
            int lower = mod * mod;
            int upper = (int)Math.Pow(mod, 1 + Degree() / 2);
            for (int n = lower; n < upper; n++)
            {
                Polynomial g = IntToPolynomial(n, mod);
                if (g.IsMonic())
                {
                    g = (Polynomial)Mod(g, mod);
                    if (g.IsZero())
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        public override string Show()
        {
            string s = "";
            for (int i = Degree(); i >= 1; i--)
            {
                if (c[i] != 0)
                {
                    if (c[i] == 1)
                    {
                        s += "X";
                    }
                    else
                    {
                        s += c[i].ToString() + "X";
                    }
                    if (i != 1)
                    {
                        s += "^" + i.ToString();
                    }
                    s += " + ";
                }
            }
            s += c[0] + ".";
            return s;
        }
    }
}
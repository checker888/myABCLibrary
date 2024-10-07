import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Deque;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Objects;
import java.util.PriorityQueue;
import java.util.Scanner;


//    int n = Integer.parseInt(sc.next());
//    int [] a = arrayInputInt(n);
//    int [][] a = arrayInputInt(y, x);
//    String s = sc.next();
//    String w [] = s.split("");
//    List<Integer> list = new ArrayList<Integer>();
//    HashMap<String,Integer> map = new HashMap<String,Integer>();
//    LinkedHashMap<String,Integer> map = new LinkedHashMap<String,Integer>();

//    Deque<Integer> dq = new ArrayDeque<>();
//    TreeSet<Integer> set = new TreeSet<Integer>();
//    PriorityQueue<Integer> pq = new PriorityQueue<Integer>(Collections.reverseOrder());

//    setAdList(n);// adlist(隣接リスト)を生成
//    setPairList(n);// pairlist(ペア隣接リスト)を生成(ダイクストラ法使うとき)

//    for(String key : map.keySet()){}

//    System.out.println(String.format("%.1f", 21.8755));
//    arrayPrint(a, 0);//配列出力　空白なし

//ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
public class Main {
    static Scanner sc;
    static PrintWriter pw = new PrintWriter(System.out);
    
    static int mod = 1000000007;
    static int bigint = 2000000000;
    static long biglong = 2000000000000000000L;
    static int visited[];
    static ArrayList<ArrayList<Pair<Integer,Long>>> pairlist;
    static long cur[];
    static PriorityQueue<Pair<Long,Integer>> dijkpq = new PriorityQueue<Pair<Long,Integer>>();
    public static void main(String[] args) throws Exception{//ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
//        sc = new Scanner(       new File("src/data.txt")         );
        sc = new Scanner(       System.in       );
//        sc = new FastScanner();
        
    }//ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー

    
    
    
    //グラフーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
    
    public static void dfs(Graph g, int pos) {
        g.visited[pos] = 1;
        int size = g.adlist.get(pos).size();
        for(int i=0;i<size;i++) {
            int to = g.adlist.get(pos).get(i);
            if(g.visited[to] == 0) {
                
                dfs(g, to);
                
            }
        }
    }
    
    
    public static void dfsStack(Graph g,int start) {
        if(g.visited[start] == 0) {
            g.visited[start] = 1;
            g.dfsst.push(start);
        }
        while(!g.dfsst.isEmpty()) {

            int pos = g.dfsst.pop();
            int size = g.adlist.get(pos).size();
            for(int i=0;i<size;i++) {
                int to = g.adlist.get(pos).get(i);
                if(g.visited[to] == 0) {
                    g.visited[to] = 1;
                    //ここに処理入れる
                    g.dfsst.push(to);
                    
                }
            }
        }

    }
    
    
    //幅優先探索
    public static void bfs(Graph g,int pos) {
        
        int nn = g.visited.length;
        for(int i=0;i<nn;i++) g.visited[i] = -1;
        
        g.bfsq.add(pos);
        g.visited[pos] = 0;
        while(!g.bfsq.isEmpty()) {
            pos = g.bfsq.poll();
            int size = g.adlist.get(pos).size();
            for(int i=0;i<size;i++) {
                int to = g.adlist.get(pos).get(i);
                if(g.visited[to] == -1) {
                    g.visited[to] = g.visited[pos] +1;
                    g.bfsq.add(to);
                }
            }
        }
    }
    
    public static void setPairList(int n) {//長さが負の辺があるとき、ダイクストラ法は正しく動かない
        pairlist = new ArrayList<ArrayList<Pair<Integer,Long>>>();
        for(int i=0; i<=n; i++){
            ArrayList<Pair<Integer,Long>> ar = new ArrayList<Pair<Integer,Long>>();
            pairlist.add(ar);   
            
        }
        visited = new int[n+1];
        cur = new long[n+1];
        for(int i=0;i<=n;i++) cur[i] = biglong;
    }
    
    public static void dijkstra(int pos) {//最初の頂点posから各頂点までの距離を求める
        cur[pos] =0;
        dijkpq.add(new Pair<Long, Integer>(cur[pos], pos));
        
        while(!dijkpq.isEmpty()) {
            pos = dijkpq.poll().getRight();//優先キュー先頭の頂点をとりだす
            if(visited[pos] == 1) continue;
            
            visited[pos] = 1;
            for(int i = 0;i < pairlist.get(pos).size();i++) {
                int to = pairlist.get(pos).get(i).getLeft();
                long cost = pairlist.get(pos).get(i).getRight();
                if(cur[to] > cur[pos]+cost) {
                    cur[to] = cur[pos]+cost;
                    dijkpq.add(new Pair<Long, Integer>(cur[to], to));
                }
            }
        }
    }
    //探索ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
    //順列全探索ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー

    public static void procPerm(int[] perm) {
        //ここに処理いれる
//        System.out.println(Arrays.toString(perm));

        
        
    }
    
    public static void permutation(int[] seed) {
        int[] perm = new int[seed.length];
        boolean[] used = new boolean[seed.length];
        buildPerm(seed, perm, used, 0);
    }

    public static void buildPerm(int[] seed, int[] perm, boolean[] used, int index) {
        if (index == seed.length) {
            procPerm(perm);
            return;
        }

        for (int i = 0; i < seed.length; i++) {
            if (used[i])
                continue;
            perm[index] = seed[i];
            used[i] = true;
            buildPerm(seed, perm, used, index + 1);
            used[i] = false;
        }
    }

    
    //bit全探索のやり方思い出す用
    public static void bitSearch(int n) {
        for (int bit = 0; bit < (1 << n); bit++) {
            
            for (int i = 0; i < n; i++) {
                System.out.println(bit +" "+ (1 << i)+" "+ (bit & (1 << i)));
                
                if((bit & (1 << i)) != 0) {
                    //ここに処理いれる
                }

            }

        }
    }
    
    //指定した値以上の値が最初に出てくるidx
    public static int lowerBound(int []input, int value) {
        int low = -1;
        int high = input.length;
        
        while(high - low > 1) {
            int mid = (low + high)/2;
            if(input[mid] >= value) high = mid;
            else low = mid;
        }
        return high;
    }
    //指定した値より大きい値が最初に出てくるidx
    public static int upperBound(int []input, int value) {
        value++;
        int low = -1;
        int high = input.length;
        
        while(high - low > 1) {
            int mid = (low + high)/2;
            if(input[mid] >= value) high = mid;
            else low = mid;
        }
        return high;
    }
    
    //入出力系---------------------------------------------------------------------------------------------
    //1次元配列の入力---------------------------------------------------------------------------------------------
    public static int [] arrayInputInt(int n) {
        int [] a = new int[n];
        
        Arrays.setAll(a, i -> Integer.parseInt(sc.next()));
        return a;
    }
    
    public static long [] arrayInputLong(int n) {
        long [] a = new long[n];
        Arrays.setAll(a, i -> sc.nextLong());
        return a;
    }
    
    public static double [] arrayInputDouble(int n) {
        double [] a = new double[n];
        Arrays.setAll(a, i -> Double.parseDouble(sc.next()));
        return a;
    }
    
    public static String [] arrayInputString(int n) {
        String [] a = new String[n];
        Arrays.setAll(a, i -> sc.next());
        return a;
    }
    
    
    //1次元累積和の作成
    public static int [] csum(int []input) {
        int n = input.length;
        int [] csum = new int[n+1];
        csum[0] = 0;
        for(int i = 1;i < n+1;i++){
            csum[i] = csum[i-1] + input[i-1];
        }
        return csum;
    }
    public static long [] csum(long []input) {
        int n = input.length;
        long [] csum = new long[n+1];
        csum[0] = 0;
        for(int i = 1;i < n+1;i++){
            csum[i] = csum[i-1] + input[i-1];
        }
        return csum;
    }
    
    
    //2次元配列の入力---------------------------------------------------------------------------------------------
    public static int [][] arrayInputInt(int y, int x) {
        int [][] a = new int[y][x];
        for(int i = 0;i < y;i++){
            for(int j = 0;j < x;j++) {
                a[i][j] = Integer.parseInt(sc.next());
            }
        }
        return a;
    }
    
    public static long [][] arrayInputLong(int y, int x) {
        long [][] a = new long[y][x];
        for(int i = 0;i < y;i++){
            for(int j = 0;j < x;j++) {
                a[i][j] = Long.parseLong(sc.next());
            }
        }
        return a;
    }
    
    public static double [][] arrayInputDouble(int y, int x) {
        double [][] a = new double[y][x];
        for(int i = 0;i < y;i++){
            for(int j = 0;j < x;j++) {
                a[i][j] = Double.parseDouble(sc.next());
            }
        }
        return a;
    }
    
    public static String [][] arrayInputString(int y, int x) {
        String [][] a = new String[y][x];
        for(int i = 0;i < y;i++){
            for(int j = 0;j < x;j++) {
                a[i][j] = sc.next();
            }
        }
        return a;
    }
    
    public static int [][] csum(int [][]input) {
        int h = input.length;
        int w = input[0].length;
        int [][] csum = new int[h+1][w+1];
        for(int i = 0;i < h+1;i++){
            for(int j = 0;j< w+1;j++) {
                if(i==0 || j==0) csum[i][j] = 0;
                else csum[i][j] = csum[i-1][j] + csum[i][j-1] - csum[i-1][j-1] + input[i-1][j-1];
            }
            
        }
        return csum;
    }
    
    public static long [][] csum(long [][]input) {
        int h = input.length;
        int w = input[0].length;
        long [][] csum = new long[h+1][w+1];
        for(int i = 0;i < h+1;i++){
            for(int j = 0;j< w+1;j++) {
                if(i==0 || j==0) csum[i][j] = 0;
                else csum[i][j] = csum[i-1][j] + csum[i][j-1] - csum[i-1][j-1] + input[i-1][j-1];
            }
            
        }
        return csum;
    }
    
    //1次元配列の出力(blank  0:間の空白なし　1:あり)---------------------------------------------------------------------------------------------

    public static void arrayPrint(int[] n, int blank) {
        for(int x = 0;x < n.length; x++) {
            
            pw.print(n[x]);
            if(x != n.length -1 && blank == 1) pw.print(" ");
            
        }
        pw.println();
        pw.flush();
    }
    
    public static void arrayPrint(long[] n, int blank) {
        for(int x = 0;x < n.length; x++) {
            
            pw.print(n[x]);
            if(x != n.length -1 && blank == 1) pw.print(" ");
            
        }
        pw.println();
        pw.flush();
    }
    
    public static void arrayPrint(double[] n, int blank) {
        for(int x = 0;x < n.length; x++) {
            
            pw.print(n[x]);
            if(x != n.length -1 && blank == 1) pw.print(" ");
            
        }
        pw.println();
        pw.flush();
    }
    
    public static void arrayPrint(String[] n, int blank) {
        for(int x = 0;x < n.length; x++) {
            
            pw.print(n[x]);
            if(x != n.length -1 && blank == 1) pw.print(" ");
            
        }
        pw.println();
        pw.flush();
    }
    
    public static void arrayPrint(boolean[] n, int blank) {
        for(int x = 0;x < n.length; x++) {
            if(n[x]) pw.print("T");
            else pw.print("F");
            
            if(x != n.length -1 && blank == 1) pw.print(" ");
            
        }
        pw.println();
        pw.flush();
    }
    
    //2次元配列の出力
    public static void arrayPrint(int[][] n, int blank) {
//        System.out.println(n.length+" "+n[1].length);
        for(int y = 0 ;y < n.length ; y++) {
            for(int x = 0;x < n[0].length; x++) {
                pw.print(n[y][x]);
                if(x != n[0].length -1 && blank == 1) pw.print(" ");
            }
            pw.println();
        }
        pw.flush();
    }
    
    public static void arrayPrint(long[][] n, int blank) {
        for(int y = 0 ;y < n.length; y++) {
            for(int x = 0;x < n[0].length; x++) {
                pw.print(n[y][x]);
                if(x != n[0].length -1 && blank == 1) pw.print(" ");
            }
            pw.println();
        }
        pw.flush();
    }
    
    public static void arrayPrint(double[][] n, int blank) {
        for(int y = 0 ;y < n.length; y++) {
            for(int x = 0;x < n[0].length; x++) {
                pw.print(n[y][x]);
                if(x != n[0].length -1 && blank == 1) pw.print(" ");
            }
            pw.println();
        }
        pw.flush();
    }
    
    public static void arrayPrint(String[][] n, int blank) {
        for(int y = 0 ;y < n.length; y++) {
            for(int x = 0;x < n[0].length; x++) {
                pw.print(n[y][x]);
                if(x != n[0].length -1 && blank == 1) pw.print(" ");
            }
            pw.println();
        }
        pw.flush();
    }
    
    public static void arrayPrint(boolean[][] n, int blank) {
      for(int y = 0 ;y < n.length ; y++) {
          for(int x = 0;x < n[0].length; x++) {
              if(n[y][x]) pw.print("T");
              else pw.print("F");
              if(x != n[0].length -1 && blank == 1) pw.print(" ");
          }
          pw.println();
      }
      pw.flush();
  }
   
    
    //ーーーーーーーーーーーーーーーーーーーーーーーーーー
    
    //数値を1桁ずつ取り出してarraylistに格納（向き逆---------------------------------------------------------------------------------------------
    public static List<Integer> cutInt(int n) {
        List<Integer> list = new ArrayList<>();
        
        while(n != 0) {
            list.add(n%10);
            n /=10;
        }
        
        return list;
    }
    
    //配列の値を大きい順に並べ替える---------------------------------------------------------------------------------------------
    public static int [] sortDesc(int []n) {
        Arrays.sort(n);
        int [] m = new int[n.length];
        for(int i=0;i<n.length;i++) {
            m[i] = n[n.length - i - 1];
        }
        return m;
    }
    
    
    
    
    //int、long系
    //10進数のintを2進数のStringで返す　値が65535より大きくなるときはintにできないので注意
    public static String dtob (int d) {//Decimal to Binary
        StringBuilder b = new StringBuilder();
        if(d == 0) return "0";
        
        while(d != 0) {
            int work = d % 2;
            b.insert(0, work);
            d /= 2;
        }
        
        return b.toString();
    }
    
    
    //3つの数のmax
    public static int max3(int a, int b, int c) {
        return Math.max(a, Math.max(b, c));
    }
    //3つの数のmin
    public static int min3(int a, int b, int c) {
        return Math.min(a, Math.min(b, c));
    }
    
    
    //double系---------------------------------------------------------------------------------------------
    
  //少数点のf桁までを出力（f+1桁を四捨五入)
    public static void doublePrint(double n, int f) {
        Integer i = Integer.valueOf(f);
        String s = "%."+i.toString()+"f";
        System.out.println(String.format(s, n));
    }
    
    
    
    
    
    
    
    
    
    
    
    //String系---------------------------------------------------------------------------------------------
    //wの中から、連続するk文字の部分文字列が回文かどうか調べる
    public static boolean isPalindrome(String w[], int k) {
        
        boolean f = false;
        int n = w.length;
        Outer:
        for(int i=0;i<=n-k;i++) {
            for(int j=1;j<=k/2;j++) {
                if(!w[i+j-1].equals(w[i+k-j])) {
                    break;
                }else if(j == k/2) {
                    f =true;
                    break Outer;
                }
            }
            
        }
        return f;
    }
    
    
    
    //hashmap系---------------------------------------------------------------------------------------------

    
    
    //アルファベットをカウントするhashmap　小文字
    public static LinkedHashMap<String,Integer> lowerABMap () {
        LinkedHashMap<String,Integer> map = new LinkedHashMap<String,Integer>();
        for(int i=0;i<26;i++) {
            char c = (char)('a'+i);
//            System.out.println(c);
            map.put(String.valueOf(c), 0);
        }
        return map;
    }
    
    
    //アルファベットをカウントするhashmap　大文字
    public static LinkedHashMap<String,Integer> upperABMap () {
        LinkedHashMap<String,Integer> map = new LinkedHashMap<String,Integer>();
        for(int i=0;i<26;i++) {
            char c = (char)('A'+i);
//            System.out.println(c);
            map.put(String.valueOf(c), 0);
        }
        return map;
    }
    
    
    
    //数学系ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
    //階乗(!n)
    static long factorial(int n) {
        long ans = 1;
        for(int i =2; i <= n; i++) {
            ans *=(long)i;
        }
        return ans;
    }
    //順列(nPk)
    static long permutation(int n, int k) {
        long ans = 1;
        for(int i = n - k + 1; i <= n; i++) {
            ans *=(long)i;
        } 
        return ans;
    }
    //組合せ(nCk)
    static long combination(int n, int k) {
        return permutation(n, k) / factorial(k);
    }
    
    //素数判定( O(ルートn) )
    static boolean isPrime(int n) {
        if(n != 2) {
            if(n % 2 == 0) return false;
        }
        for(int i = 3; i * i <= n; i+=2) {
            if(n % i == 0) return false;
        }
        return true;
    }
    
    //最大公約数求める
    static long gcd(long a, long b) {
        while(a>=1 && b>=1) {
            if(a >= b) a %= b;
            else b %= a;
        }
        if(a != 0) return a;
        else return b;
    }
    //最小公倍数求める
    static long lcm(long a, long b) {
        return(a * b / gcd(a, b));
    }
    //aのb乗求める
    static long power(long a, long b) {
        long p = a;
        long ans = 1;
        int lim = 30;
        for(int i=0;i<lim;i++) {
            int wari = (1 << i);
            if((b/wari) % 2 == 1) {
                ans= (ans * p);
            }
            p = (p*p);
        }
        return ans;
    }
    //aのb乗をmodで割る
    static long power(long a, long b, long m) {
        long p = a;
        long ans = 1;
        int lim = 30;
        for(int i=0;i<lim;i++) {
            int wari = (1 << i);
            if((b/wari) % 2 == 1) {
                ans= (ans * p) % m;
            }
            p = (p*p)%m;
        }
        return ans;
    }
    
    //aの逆元求める(modで割る
    static long modinv(long a, long m) {
        long b = m;
        long u = 1;
        long v = 0;
        while(b > 0) {
            long t = a/b;
            
            a -= t*b; 
            long work = a;
            a = b;
            b = work;
            
            u -= t*v;
            work = u;
            u = v;
            v = work;
        }
        
        u %= m;
        if(u < 0) u+=m;
        return u;
    }
    
    //nの約数列挙
    public static void printDiv(long n) {
        Deque<Long> dq = new ArrayDeque<>();
        for(long i=1;i*i<=n;i++) {
            if(n%i == 0) {
                pw.println(i);
                dq.push(n/i);
            }
        }
        while(!dq.isEmpty()) {
            pw.println(dq.pop());
        }
        pw.flush();
    }
    


}

class Edge{
    long u;
    long v;
    long cost;
}


class Graph {
    ArrayList<ArrayList<Integer>> adlist = new ArrayList<ArrayList<Integer>>();//隣接リスト
    int n;
    int visited[];
    Deque<Integer> dfsst = new ArrayDeque<>();//dfs用スタック
    Deque<Integer> bfsq = new ArrayDeque<>();//bfs用キュー
    UnionFind uf;
    //コンストラクタ
    public Graph(int n) {
        this.n = n;
        for(int i=0; i<=n; i++){
            ArrayList<Integer> ar = new ArrayList<Integer>();   //外側のList生成
            this.adlist.add(ar);   
            
        }
        this.visited = new int[n+1];
        uf = new UnionFind(n+1);
    }
    
    //単一方向（uからv）に頂点を接続
    public void connect(int u, int v) {
        this.adlist.get(u).add(v);
        uf.unite(u, v);
        }
    //双方向に頂点を接続
    public void connectMutual(int u, int v){
        this.adlist.get(u).add(v);
        this.adlist.get(v).add(u);
        uf.unite(u, v);
    }
    //unionfindで2つの頂点が繋がっているか確認
    public boolean isConnect(int u, int v) {return uf.same(u, v);}
    
    //探索済み（または最短距離）を示すvisitedのリセット
    public void resetVisited() {this.visited = new int[n+1];}


}


class UnionFind {
    int[] parent;
    int[] rank;
    //コンストラクタ
    public UnionFind(int n) {
        this.parent = new int[n];
        this.rank = new int[n];
        // 最初はすべてが根(or根はない)
        for (int i = 0; i < n; i++) {
            parent[i] = i;
//            parent[i] = -1; //根はない
            rank[i] = 0;
        }
    }
//     要素の根を経路圧縮して返す
    public int root(int x) {
        if(x != parent[x]) {
            parent[x] = root(parent[x]);
        }
        return parent[x];
    }
    //２つの要素が同じ集合に属するかどうかを返す
    public boolean same(int x, int y) {
        return root(x) == root(y);
    }
    //xの属する集合とyの属する集合をつなげる
    public void unite(int x, int y) {
        x = root(x);
        y = root(y);

        if (x == y)  return;// 属する集合が同じな場合、何もしない

        if (rank[x] > rank[y]) {
            parent[y] = x;
        } else {
            parent[x] = y;
            if(rank[x] == rank[y]) rank[y]++;
        } 
    }
}

class Pair<S extends Comparable<S>, T extends Comparable<T>> implements Comparable<Pair<S,T>>{//ソートできるようにComparable使う
    S left;
    T right;

    public Pair(S s, T t){
        left = s;
        right = t;
    }
    //getメソッドは、必ず代入してから使う　（if文などで比較したとき、intをInteger型で比較して失敗したりすることに気づかない）
    public S getLeft(){return left;}
    public T getRight(){return right;}

    public int compareTo(Pair<S,T> another){
        java.util.Comparator<Pair<S,T>> comp1 = java.util.Comparator.comparing(Pair::getLeft);
        java.util.Comparator<Pair<S,T>> comp2 = comp1.thenComparing(Pair::getRight);
        return comp2.compare(this, another);
    }

}

import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
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

//    Deque<Integer> dq = new ArrayDeque<>();
//    TreeSet<Integer> set = new TreeSet<Integer>();
//    PriorityQueue<Integer> pq = new PriorityQueue<Integer>(Collections.reverseOrder());

//    for(String key : map.keySet()){}
//    Arrays.sort(aw,Collections.reverseOrder());
//    System.out.println(String.format("%.1f", 21.8755));

//ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
public class Main {
    static Scanner sc;
    static PrintWriter pw = new PrintWriter(System.out);
    
    static int mod = 1000000007;
    static int bigint = 2000000000;
    static long biglong = 2000000000000000000L;
    public static void main(String[] args) throws Exception{//ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
//        sc = new Scanner(       new File("src/data.txt")         );
        sc = new Scanner(       System.in       );
//        sc = new FastScanner();
        
    }//ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー

    
    
    
    //グラフーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
    
    //深さ優先探索
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
        Deque<Integer> st = new ArrayDeque<>();//dfs用スタック
        if(g.visited[start] == 0) {
            g.visited[start] = 1;
            st.push(start);
        }
        while(!st.isEmpty()) {

            int pos = st.pop();
            int size = g.adlist.get(pos).size();
            for(int i=0;i<size;i++) {
                int to = g.adlist.get(pos).get(i);
                if(g.visited[to] == 0) {
                    g.visited[to] = 1;
                    //ここに処理入れる
                    st.push(to);
                    
                }
            }
        }

    }
    
    
    //幅優先探索
    public static void bfs(Graph g,int pos) {
        Deque<Integer> bfsq = new ArrayDeque<>();//bfs用キュー
        int nn = g.visited.length;
        for(int i=0;i<nn;i++) g.visited[i] = -1;//計算量注意！！！！！！！！！！！！
        
        bfsq.add(pos);
        g.visited[pos] = 0;
        while(!bfsq.isEmpty()) {
            pos = bfsq.poll();
            int size = g.adlist.get(pos).size();
            for(int i=0;i<size;i++) {
                int to = g.adlist.get(pos).get(i);
                if(g.visited[to] == -1) {
                    g.visited[to] = g.visited[pos] +1;
                    bfsq.add(to);
                }
            }
        }
    }
    //始点が複数あるときに使う多視点BFS
    public static void multi_source_bfs(Graph g,Deque<Integer> que) {
        Deque<Integer> bfsq = new ArrayDeque<>();//bfs用キュー
        int nn = g.visited.length;
        for(int i=0;i<nn;i++) g.visited[i] = -1;//計算量注意！！！！！！！！！！！！
        while(!que.isEmpty()) {
            int pos = que.poll();
            bfsq.add(pos);
            g.visited[pos] = 0;
        }

        while(!bfsq.isEmpty()) {
            int pos = bfsq.poll();
            int size = g.adlist.get(pos).size();
            for(int i=0;i<size;i++) {
                int to = g.adlist.get(pos).get(i);
                if(g.visited[to] == -1) {
                    g.visited[to] = g.visited[pos] +1;
                    bfsq.add(to);
                }
            }
        }
    }
    //ダイクストラ法
    
    public static void dijkstra(WeightedGraph g,int pos) {//最初の頂点posから各頂点までの距離を求める
        g.resetCur();
        PriorityQueue<Pair<Long,Integer>> dijkpq = new PriorityQueue<Pair<Long,Integer>>();//pair(最小コスト候補、頂点)
        g.cur[pos] =0;
        dijkpq.add(new Pair<Long, Integer>(g.cur[pos], pos));
        
        while(!dijkpq.isEmpty()) {
            pos = dijkpq.poll().getRight();//優先キュー先頭の頂点をとりだす
            if(g.visited[pos] == 1) continue;
            
            g.visited[pos] = 1;
            for(int i = 0;i < g.pairadlist.get(pos).size();i++) {
                int to = g.pairadlist.get(pos).get(i).to;
                long cost = g.pairadlist.get(pos).get(i).cost;
                if(g.cur[to] > g.cur[pos]+cost) {
                    g.cur[to] = g.cur[pos]+cost;
                    dijkpq.add(new Pair<Long, Integer>(g.cur[to], to));
                }
            }
        }
    }
    //ワーシャルフロイド法(dpを返す)
    public static long[][] WarshallFloyd(WeightedGraph g){
        long dp[][] = new long[g.n+1][g.n+1];
        //dp初期化
        for(int i=0;i<=g.n;i++) {
            for(int j=0;j<=g.n;j++) {
                if(i == j) dp[i][j] = 0;
                else dp[i][j] = biglong/2;
            }
        }
        //初期状態の設定(隣接する頂点のみのコスト)
        for(int from=0;from<=g.n;from++) {
            ArrayList<Edge> a = g.pairadlist.get(from);
            for(int j=0;j<a.size();j++) {
                Edge p = a.get(j);
                int to = p.to;
                dp[from][to] = p.cost;
                
            }
        }
//        arrayPrint(dp,1);
        //ワーシャルフロイド法
        for(int k=0;k<=g.n;k++){
            for(int i=0;i<=g.n;i++){
                for(int j=0;j<=g.n;j++)  dp[i][j]=Math.min(dp[i][j],dp[i][k]+dp[k][j]);
            }
        }
        return dp;
    }
    
    //クラスカル法
    //最小全域木の構造ごと返す
    public static WeightedGraph kruskal(WeightedGraph g) {
        WeightedGraph result = new WeightedGraph(g.n);
        UnionFind uf = new UnionFind(g.n+1);
        Collections.sort(g.edges, new Comparator<int[]>() {//コストの小さい順にソート
            @Override
            public int compare(int []edge1, int []edge2) {
                return edge1[2] - edge2[2];
            }
        });
        result.costsum = 0;
        for(int [] edge : g.edges) {
            int u = edge[0];
            int v = edge[1];
            long cost = edge[2];
            if(!uf.same(u, v)) {
//                result.connect(u, v, cost);//有向グラフのとき
                result.connectMutual(u, v, cost);//無向グラフのとき (コスト総和の計算も関数内で行われてる
                uf.unite(u, v);
            }
        }
        return result;
    }
    //最小全域木のコストの総和のみを記録する
//    public static void kruskal(WeightedGraph g) {
//        Collections.sort(g.edges, new Comparator<int[]>() {//コストの小さい順にソート
//            @Override
//            public int compare(int []edge1, int []edge2) {
//                return edge1[2] - edge2[2];
//            }
//        });
//        UnionFind ufk = new UnionFind(g.n+1);
//        g.costsum = 0;
//        for(int [] edge : g.edges) {
//            int u = edge[0];
//            int v = edge[1];
//            long cost = edge[2];
//            if(!ufk.same(u, v)) {
//                ufk.unite(u, v);
//                g.costsum += cost;
//            }
//        }
//    }
    

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
            if (used[i])continue;
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
//                System.out.println(bit +" "+ (1 << i)+" "+ (bit & (1 << i)));
                if((bit & (1 << i)) != 0) {
                    
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
    //valueの個数取得
    public static int count(int[] input, int value) {
        return upperBound(input, value) - lowerBound(input, value);
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
    public static void arrayPrint(int[] n) {arrayPrint(n,1);}
    public static void arrayPrint(long[] n) {arrayPrint(n,1);}
    public static void arrayPrint(double[] n) {arrayPrint(n,1);}
    public static void arrayPrint(String[] n) {arrayPrint(n,0);}
    public static void arrayPrint(boolean[] n) {arrayPrint(n,1);}
    public static void arrayPrint(Object[] n) {arrayPrint(n,1);}
    public static void arrayPrint(int[][] n) {arrayPrint(n,1);}
    public static void arrayPrint(long[][] n) {arrayPrint(n,1);}
    public static void arrayPrint(double[][] n) {arrayPrint(n,1);}
    public static void arrayPrint(String[][] n) {arrayPrint(n,0);}
    public static void arrayPrint(boolean[][] n) {arrayPrint(n,1);}
    public static void arrayPrint(Object[][] n) {arrayPrint(n,1);}
    
    
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
    
    
    public static void arrayPrint(Object[] n, int blank) {
        for(int x = 0;x < n.length; x++) {

            pw.print(n[x]);
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
    public static int max3(int a, int b, int c) {return Math.max(a, Math.max(b, c));}
    public static int min3(int a, int b, int c) {return Math.min(a, Math.min(b, c));}
    public static int max4(int a, int b, int c,int d) {return Math.max(a, Math.max(b, Math.max(c, d)));}
    public static int min4(int a, int b, int c,int d) {return Math.min(a, Math.min(b, Math.min(c, d)));}
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
        while (b != 0) {
            long tmp = a % b;
            a = b;
            b = tmp;
        }
        return a;
    }
    //最小公倍数求める
    static long lcm(long a, long b) {
        return(a * b / gcd(a, b));
    }
    //正確な平方根をとる
    public static long sqrt(long x) {
        long r = (long)Math.sqrt(x);
        while (r * r > x) r--;
        while ((r + 1) * (r + 1) <= x) r++;
        return (r * r == x) ? r : -1;
    }
    //aのb乗求める
    static long power(long a, long b) {
        long ans = 1;
        while (b > 0) {
            if ((b & 1) == 1) ans *= a; //ビットが立っているときだけ掛ける
            a *= a;                     //底を二乗
            b >>= 1;                    //指数を1bit右シフト
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
    
    public static boolean[] sieveOfEratosthenes(int n) {
        boolean[] isPrime = new boolean[n+1];
        for(int i=0;i<=n;i++)isPrime[i] = true;
        isPrime[0] = isPrime[1] = false;
        
        for(int p=2;p<=n;p++) {
            if(!isPrime[p]) continue;
            
            for(int q = p*2;q<=n;q+=p) {
                isPrime[q] = false;
            }
            
        }
        return isPrime;
    }
    
    public static double log2(long x) { return Math.log(x) / Math.log(2); }
    public static double log10(long x) { return Math.log10(x); }
    
    //10進数のdをp進数のStringで返す　
    public static String decimalConversion (long d, int p) {
        StringBuilder b = new StringBuilder();
        if(d == 0) return "0";

        while(d != 0) {
            long work = d % p;
            b.insert(0, work);
            d /= p;
        }

        return b.toString();
    }
  //from進数の数であるnをto進数に変換
    public static String baseConversion(String n, int from, int to) {
        String nw[] = n.split("");
        long p = 1;
        long sum10 = 0;
        //10進数に変換
        for(int i=nw.length-1; i>=0; i--) {
            long a = Long.valueOf(nw[i]) * p;
            sum10 += a;
            p *= from;
        }

        String cnvSum = decimalConversion(sum10, to);
        return cnvSum;
    }
}

//重みなしグラフ
class Graph {
  ArrayList<ArrayList<Integer>> adlist = new ArrayList<ArrayList<Integer>>();//隣接リスト
  int n;
  int visited[];
  //コンストラクタ
  public Graph(int n) {
      this.n = n;
      for(int i = 0; i <= n; i++) adlist.add(new ArrayList<>());
      this.visited = new int[n+1];
  }
  
  //単一方向（uからv）に頂点を接続
  public void connect(int u, int v) {
      this.adlist.get(u).add(v);
//      if (!adlist.get(u).contains(v)) adlist.get(u).add(v);
  }
  //双方向に頂点を接続
  public void connectMutual(int u, int v){
      this.adlist.get(u).add(v);
      this.adlist.get(v).add(u);
  }
  //unionfindで2つの頂点が繋がっているか確認
//  public boolean isConnect(int u, int v) {return this.uf.same(u, v);}
  //探索済み（または最短距離）を示すvisitedのリセット
  public void resetVisited(int value) { Arrays.fill(this.visited, value);}
}

//重みつきグラフ
class WeightedGraph {
  static long biglong = 2000000000000000000L;
  ArrayList<ArrayList<Edge>> pairadlist;//pair(移動先の頂点to、コストcost)
  
  List<int[]> edges ;
  int n;
  int visited[];
  long cur[];//ダイクストラ法で使う、頂点までのコストの総和の最小値
  long costsum;//クラスカル法で使う、すべての辺のコストの総和
  //コンストラクタ
  public WeightedGraph(int n) {
      pairadlist = new ArrayList<ArrayList<Edge>>();
//      edges = new ArrayList<int[]>();//クラスカル法のときのみ
      this.n = n;
      for(int i=0; i<=n; i++)this.pairadlist.add(new ArrayList<Edge>());
      this.visited = new int[n+1];
      this.cur = new long[n+1];
  }
  
  //単一方向（uからv）に頂点を接続
  public void connect(int u, int v, long cost) {
      this.pairadlist.get(u).add(new Edge(v,cost));
      this.costsum += cost;
//      edges.add(new int[] { u, v ,(int)cost });//クラスカル法のときのみ
  }
  
  //双方向に頂点を接続
  public void connectMutual(int u, int v, long cost){
      this.pairadlist.get(u).add(new Edge(v, cost));
      this.pairadlist.get(v).add(new Edge(u, cost));
      this.costsum += cost;
//      edges.add(new int[] { u, v ,(int)cost });//クラスカル法のときのみ
  }
  //unionfindで2つの頂点が繋がっているか確認
//  public boolean isConnect(int u, int v) {return this.uf.same(u, v);}
  //探索済み（または最短距離）を示すvisitedのリセット
  public void resetVisited(int value) { Arrays.fill(this.visited, value);}
  
  public void resetCur() {Arrays.fill(cur, biglong);}
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
//          parent[i] = -1; //根はない
          rank[i] = 0;
      }
  }
//   要素の根を経路圧縮して返す
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

class Edge {
  int to;
  long cost;
  Edge(int to, long cost) {
      this.to = to;
      this.cost = cost;
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
    
    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }

        @SuppressWarnings("unchecked")
        Pair<S, T> another = (Pair<S, T>) obj;
        return Objects.equals(this.getLeft(), another.getLeft())
                && Objects.equals(this.getRight(), another.getRight());

    }

    @Override
    public int hashCode() {
        int h1 = left.hashCode();
        int h2 = right.hashCode();
        String connected = String.valueOf(h1) + ' ' + String.valueOf(h2);
        return connected.hashCode();
    }
    @Override
    public String toString() {
        return left+" "+right;
    }
}


class Triple<S extends Comparable<S>, T extends Comparable<T>, U extends Comparable<U>> implements Comparable<Triple<S,T,U>>{//ソートできるようにComparable使う
    S left;
    T mid;
    U right;
    public Triple(S s, T t, U u){
        left = s;
        mid = t;
        right = u;
    }
    //getメソッドは、必ず代入してから使う　（if文などで比較したとき、intをInteger型で比較して失敗したりすることに気づかない）
    public S getLeft(){return left;}
    public T getMid(){return mid;}
    public U getRight() {return right;}
    
    public int compareTo(Triple<S,T,U> another){
        java.util.Comparator<Triple<S,T,U>> comp1 = java.util.Comparator.comparing(Triple::getLeft);
        java.util.Comparator<Triple<S,T,U>> comp2 = comp1.thenComparing(Triple::getMid);
        java.util.Comparator<Triple<S,T,U>> comp3 = comp2.thenComparing(Triple::getRight);
        return comp3.compare(this, another);
    }
    
    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }

        @SuppressWarnings("unchecked")
        Triple<S, T, U> another = (Triple<S, T, U>) obj;
        return Objects.equals(this.getLeft(), another.getLeft())
                && Objects.equals(this.getRight(), another.getRight())
                    && Objects.equals(this.getMid(), another.getMid());

    }

    @Override
    public int hashCode() {
        int h1 = left.hashCode();
        int h2 = right.hashCode();
        int h3 = mid.hashCode();
        String connected = String.valueOf(h1) + ' ' + String.valueOf(h2) + ' ' + String.valueOf(h3);
        return connected.hashCode();
    }
    
    @Override
    public String toString() {
        return left+" "+mid+" "+right;
    }

}

class Four<S extends Comparable<S>, T extends Comparable<T>, U extends Comparable<U>, V extends Comparable<V>> implements Comparable<Four<S,T,U,V>>{//ソートできるようにComparable使う
    S left;
    T lmid;
    U rmid;
    V right;
    public Four(S s, T t, U u,V v){
        left = s;
        lmid = t;
        rmid = u;
        right = v;
    }
    //getメソッドは、必ず代入してから使う　（if文などで比較したとき、intをInteger型で比較して失敗したりすることに気づかない）
    public S getLeft(){return left;}
    public T getLMid(){return lmid;}
    public U getRMid() {return rmid;}
    public V getRight() {return right;}
    
    public int compareTo(Four<S,T,U,V> another){
        java.util.Comparator<Four<S,T,U,V>> comp1 = java.util.Comparator.comparing(Four::getLeft);
        java.util.Comparator<Four<S,T,U,V>> comp2 = comp1.thenComparing(Four::getLMid);
        java.util.Comparator<Four<S,T,U,V>> comp3 = comp2.thenComparing(Four::getRMid);
        java.util.Comparator<Four<S,T,U,V>> comp4 = comp3.thenComparing(Four::getRight);
        return comp4.compare(this, another);
    }
    
    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }

        @SuppressWarnings("unchecked")
        Four<S, T, U, V> another = (Four<S, T, U, V>) obj;
        return Objects.equals(this.getLeft(), another.getLeft())
                && Objects.equals(this.getRight(), another.getRight())
                    && Objects.equals(this.getLMid(), another.getLMid())
                        && Objects.equals(this.getRMid(), another.getRMid());

    }

    @Override
    public int hashCode() {
        int h1 = left.hashCode();
        int h2 = right.hashCode();
        int h3 = lmid.hashCode();
        int h4 = rmid.hashCode();
        String connected = String.valueOf(h1) + ' ' + String.valueOf(h2) + ' ' + String.valueOf(h3) + ' ' + String.valueOf(h4);
        return connected.hashCode();
    }
    
    @Override
    public String toString() {
        return left+" "+lmid+" "+rmid+" "+right;
    }
}





//配列の最大・最小値取得
//Arrays.stream([配列名]).min().getAsInt()
//Arrays.stream([配列名]).max().getAsInt()



//grid 上下左右・斜めの接続
//for(int i=0;i<h;i++) {
//    for(int j=0;j<w;j++) {
//        int num = (i*w)+j;
//        if(grid[i][j].equals("#")) {
//            if(i != 0 && grid[i-1][j].equals("#")) g.connect(num, num-w);
//            if(i != h-1 && grid[i+1][j].equals("#"))g.connect(num, num+w);
//            if(j != 0 && grid[i][j-1].equals("#"))g.connect(num, num-1);
//            if(j != w-1 && grid[i][j+1].equals("#"))g.connect(num, num+1);
//        
//            if(i != 0 && j != 0 && grid[i-1][j-1].equals("#")) g.connect(num, num-w-1);
//            if(i != 0 && j != w-1 && grid[i-1][j+1].equals("#")) g.connect(num, num-w+1);
//            if(i != h-1 && j != 0 && grid[i+1][j-1].equals("#")) g.connect(num, num+w-1);
//            if(i != h-1 && j != w-1 && grid[i+1][j+1].equals("#")) g.connect(num, num+w+1);
//         }
//    }
//}
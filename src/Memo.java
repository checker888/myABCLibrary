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
//    int bitCount = Integer.bitCount(num);
//    Arrays.sort(bx, (a, b) -> Long.compare(Math.abs(a), Math.abs(b)));






//配列の最大・最小値取得
//    Arrays.stream([配列名]).min().getAsInt()
//    Arrays.stream([配列名]).max().getAsInt()



// grid 上下左右・斜めの接続
//    for(int i=0;i<h;i++) {
//        for(int j=0;j<w;j++) {
//            int num = (i*w)+j;
//            if(grid[i][j].equals("#")) {
//                if(i != 0 && grid[i-1][j].equals("#")) g.connect(num, num-w);
//                if(i != h-1 && grid[i+1][j].equals("#"))g.connect(num, num+w);
//                if(j != 0 && grid[i][j-1].equals("#"))g.connect(num, num-1);
//                if(j != w-1 && grid[i][j+1].equals("#"))g.connect(num, num+1);
//        
//                if(i != 0 && j != 0 && grid[i-1][j-1].equals("#")) g.connect(num, num-w-1);
//                if(i != 0 && j != w-1 && grid[i-1][j+1].equals("#")) g.connect(num, num-w+1);
//                if(i != h-1 && j != 0 && grid[i+1][j-1].equals("#")) g.connect(num, num+w-1);
//                if(i != h-1 && j != w-1 && grid[i+1][j+1].equals("#")) g.connect(num, num+w+1);
//            }
//        }
//    }
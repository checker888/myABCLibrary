import java.io.PrintWriter;
import java.util.Random;
import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
public class MakeTestcase {
    static FileWriter file;
    static PrintWriter pw;
    
    static Random random = new Random();
    static String uppercharacters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    static String lowercharacters = "abcdefghijklmnopqrstuvwxyz";
    static String characters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    static String numbers = "1234567890";
    public static void main(String[] args) {
        setup();
        
        
        
        
        
        
        
        pw.close();//必須
    }
    
    
    public static void setup() {
        try {
            file = new FileWriter("src/testcase.txt");
            pw = new PrintWriter(new BufferedWriter(file));
        }catch (IOException e) {
            e.printStackTrace();
            System.exit(0);
        }
    }
    
    //文字列stの中から選んだn文字のランダムな文字列を生成
    public static String randomizeString(int n, String st) {
        StringBuilder randomStringBuilder = new StringBuilder(n);
        for (int i = 0; i < n; i++) {
            int index = random.nextInt(st.length());
            randomStringBuilder.append(st.charAt(index));
        }
        return randomStringBuilder.toString();
    }
    public static String randomizeString(int n) {return randomizeString(n,lowercharacters);}
    public static String randomizeString(String st) {return randomizeString(1,st);}
    
    //low <= x <= high の範囲の整数xを生成
    public static int randomizeInt(int low, int high) {
        return low + random.nextInt((high - low) + 1);
    }
    
    public static long randomizeLong(long low, long high) {
        return low + random.nextLong((high - low) + 1);
    }
    
    //low <= x <= high の範囲の実数xを生成(少数位:digit)
    public static double randomizeDouble(double low, double high, int digit) {
        double value = low + (high + 1 - low) * random.nextDouble(); 
        double scale = Math.pow(10, digit);                          
        return Math.floor(value * scale) / scale;                    
    }
    public static double randomizeDouble(double low, double high) {return randomizeDouble(low,high,1);}
    
    //low <= a[i] <= high の範囲の配列aを生成
    public static int[] randomizeIntArray(int n, int low, int high) {
        int a[] = new int[n];
        for(int i=0;i<n;i++) a[i] = randomizeInt(low, high);
        return a;
    }
    
    public static long[] randomizeLongArray(int n, long low, long high) {
        long a[] = new long[n];
        for(int i=0;i<n;i++) a[i] = randomizeLong(low, high);
        return a;
    }
    
    public static double[] randomizeDoubleArray(int n, double low, double high, int digit) {
        double[] a = new double[n];
        for (int i = 0; i < n; i++) {
            a[i] = randomizeDouble(low, high,digit);
        }
        return a;
    }
    public static double[] randomizeDoubleArray(int n, double low, double high) {return randomizeDoubleArray(n,low,high,1);}
            
    public static void arrayPrintln(int[] n) {for(int y = 0;y < n.length; y++) {pw.println(n[y]);}}
    public static void arrayPrintln(long[] n) {for(int y = 0;y < n.length; y++) {pw.println(n[y]);}}
    public static void arrayPrintln(double[] n) {for(int y = 0;y < n.length; y++) {pw.println(n[y]);}}
    public static void arrayPrintln(String[] n) {for(int y = 0;y < n.length; y++) {pw.println(n[y]);}}
    
    public static void arrayPrint(int[] n) {arrayPrint(n,1);}
    public static void arrayPrint(long[] n) {arrayPrint(n,1);}
    public static void arrayPrint(double[] n) {arrayPrint(n,1);}
    public static void arrayPrint(String[] n) {arrayPrint(n,0);}
    public static void arrayPrint(int[][] n) {arrayPrint(n,1);}
    public static void arrayPrint(long[][] n) {arrayPrint(n,1);}
    public static void arrayPrint(double[][] n) {arrayPrint(n,1);}
    public static void arrayPrint(String[][] n) {arrayPrint(n,0);}
    
    public static void arrayPrint(int[] n, int blank) {
        for(int x = 0;x < n.length; x++) {
            
            pw.print(n[x]);
            if(x != n.length -1 && blank == 1) pw.print(" ");
            
        }
        pw.println();
    }
    
    public static void arrayPrint(long[] n, int blank) {
        for(int x = 0;x < n.length; x++) {
            
            pw.print(n[x]);
            if(x != n.length -1 && blank == 1) pw.print(" ");
            
        }
        pw.println();
    }
    
    public static void arrayPrint(double[] n, int blank) {
        for(int x = 0;x < n.length; x++) {
            
            pw.print(n[x]);
            if(x != n.length -1 && blank == 1) pw.print(" ");
            
        }
        pw.println();
    }
    
    public static void arrayPrint(String[] n, int blank) {
        for(int x = 0;x < n.length; x++) {
            
            pw.print(n[x]);
            if(x != n.length -1 && blank == 1) pw.print(" ");
            
        }
        pw.println();
    }
    
    
    public static void arrayPrint(int[][] n, int blank) {
        for(int y = 0 ;y < n.length ; y++) {
            for(int x = 0;x < n[0].length; x++) {
                pw.print(n[y][x]);
                if(x != n[0].length -1 && blank == 1) pw.print(" ");
            }
            pw.println();
        }
    }
      
    public static void arrayPrint(long[][] n, int blank) {
        for(int y = 0 ;y < n.length; y++) {
            for(int x = 0;x < n[0].length; x++) {
                pw.print(n[y][x]);
                if(x != n[0].length -1 && blank == 1) pw.print(" ");
            }
            pw.println();
        }
    }
  
    public static void arrayPrint(double[][] n, int blank) {
        for(int y = 0 ;y < n.length; y++) {
            for(int x = 0;x < n[0].length; x++) {
                pw.print(n[y][x]);
                if(x != n[0].length -1 && blank == 1) pw.print(" ");
            }
            pw.println();
        }
    }
  
    public static void arrayPrint(String[][] n, int blank) {
        for(int y = 0 ;y < n.length; y++) {
            for(int x = 0;x < n[0].length; x++) {
                pw.print(n[y][x]);
                if(x != n[0].length -1 && blank == 1) pw.print(" ");
            }
            pw.println();
        }
    }
  
}

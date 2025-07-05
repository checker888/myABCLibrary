public class Debug {
    public static void main(String[] args) throws Exception{
        doMain();
    }
    
    
    public static void doMain() throws Exception{
        long start = System.currentTimeMillis();
        Main.main(new String[0]);
        long end = System.currentTimeMillis();
        System.out.println((end - start)  + "ms");
    }
    
    
    public static void doMain(int trials) throws Exception{
        for(int t=1;t<=trials;t++) {
            MakeTestcase.main(new String[0]);
            long start = System.currentTimeMillis();
            Main.main(new String[0]);
            long end = System.currentTimeMillis();
            System.out.println(t+"回目:"+(end - start)  + "ms");
        }

    }
}

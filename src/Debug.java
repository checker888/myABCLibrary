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
    
}

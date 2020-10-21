import de.lmu.ifi.dbs.elki.algorithm.Algorithm;
import de.lmu.ifi.dbs.elki.algorithm.clustering.subspace.DOC;
import de.lmu.ifi.dbs.elki.data.Cluster;
import de.lmu.ifi.dbs.elki.data.Clustering;
import de.lmu.ifi.dbs.elki.data.NumberVector;
import de.lmu.ifi.dbs.elki.data.model.SubspaceModel;
import de.lmu.ifi.dbs.elki.data.type.TypeUtil;
import de.lmu.ifi.dbs.elki.database.Database;
import de.lmu.ifi.dbs.elki.database.StaticArrayDatabase;
import de.lmu.ifi.dbs.elki.database.relation.Relation;
import de.lmu.ifi.dbs.elki.datasource.ArrayAdapterDatabaseConnection;
import de.lmu.ifi.dbs.elki.datasource.DatabaseConnection;
import de.lmu.ifi.dbs.elki.utilities.ClassGenericsUtil;
import de.lmu.ifi.dbs.elki.utilities.optionhandling.parameterization.ListParameterization;

import java.util.ArrayList;
import java.util.List;

public class ClusterApplier {
    public static void main(String[] args) {
        //Array for Testing purposes
        double[][] arr = {{1,1,1,1},{1,1,1,3},{9,9,9,5},{1,1,1,7},{1,1,1,9},{9,9,9,11}};
        //Cluster[] res = doc(arr,0.2D,0.8D,0.05D);
        //System.out.println(res.length);

    }

    /*
    *
        public static final OptionID ALPHA_ID = new OptionID("doc.alpha", "Minimum relative density for a set of points to be considered a cluster (|C|>=doc.alpha*|S|).");
        public static final OptionID BETA_ID = new OptionID("doc.beta", "Preference of cluster size versus number of relevant dimensions (higher value means higher priority on larger clusters).");
        public static final OptionID W_ID = new OptionID("doc.w", "Maximum extent of scattering of points along a single attribute for the attribute to be considered relevant.");
        public static final OptionID RANDOM_ID = new OptionID("doc.random-seed", "Random seed, for reproducible experiments.");
    *
    * */
    public static Cluster[] doc(double[][] data,
                                double alpha,
                                double beta,
                                double w) {
        DatabaseConnection dbc = new ArrayAdapterDatabaseConnection(data);
        // Create a database (which may contain multiple relations!)
        Database db = new StaticArrayDatabase(dbc, null);
        // Load the data into the database (do NOT forget to initialize...)
        db.initialize();


        Relation<NumberVector> vectors = db.getRelation(TypeUtil.NUMBER_VECTOR_FIELD);


        ListParameterization params2 = new ListParameterization();
        params2.addParameter(DOC.Parameterizer.ALPHA_ID, alpha);    //0.2D
        params2.addParameter(DOC.Parameterizer.BETA_ID, beta); //0.8D
        params2.addParameter(DOC.Parameterizer.W_ID, w);    //0.05D
        params2.addParameter(DOC.Parameterizer.RANDOM_ID, 15012021);

        Algorithm alg = ClassGenericsUtil.parameterizeOrAbort(DOC.class, params2);
        Clustering<SubspaceModel> clust = ((DOC) alg).run(db,vectors);

        assert clust != null;
        List<Cluster<SubspaceModel>> result = clust.getAllClusters();
        Cluster[] out = result.toArray(new Cluster[0]);
        return out;
    }

    public static double[][] getIds (Cluster<SubspaceModel>[] clusters){
        double[][] output= new double[clusters.length][0];
        for (int i = 0; i<clusters.length;i++){
            double[] ids = new double[clusters[i].getIDs().size()];
            String[] numSplit = clusters[i].getIDs().toString().replace(" ","").replace("[","").replace("]","").split(",");
            for(int j=0;j<numSplit.length;j++) {
                ids[j] = Double.parseDouble(numSplit[j]);
            }
            output[i] = ids;
        }
        return output;
    }

    public static double[][] getDims (Cluster<SubspaceModel>[] clusters){
        double[][] output= new double[clusters.length][0];
        for (int i = 0; i<clusters.length;i++){
            double[] dims = dimensionsFromString(clusters[i].getModel().getSubspace().dimensionsToString());
            output[i] = dims;
        }
        return output;
    }

    public static double[][] getAverages (Cluster<SubspaceModel>[] clusters){
        double[][] output= new double[clusters.length][0];
        for (int i = 0; i<clusters.length;i++){
            double[] averageValuesInCluster = clusters[i].getModel().getPrototype();
            output[i] = averageValuesInCluster;
        }
        return output;
    }
	
	public static String[] getPrototypeTypes (Cluster<SubspaceModel>[] clusters){
        String[] output= new String[clusters.length];
        for (int i = 0; i<clusters.length;i++){
            output[i] = clusters[i].getModel().getPrototypeType();
        }
        return output;
    }

    private static double[] arrayToDouble(ArrayList<Integer> ints){
        double[] result = new double[ints.size()];
        for(int i = 0;i<ints.size();i++){
            result[i]=ints.get(i);
        }
        return result;
    }

    private static double[] dimensionsFromString(String dims){
        String[] numSplit = dims.replace("[","").replace("]","").replace(" ","").split(",");
        double[] result = new double[numSplit.length];
        for(int i = 0;i<numSplit.length;i++){
            result[i]= Double.parseDouble(numSplit[i]);
        }
        return result;
    }

    public static double[][] matrix_from_array(double[] arr,int ncol) {
        if(arr.length % ncol!=0) {
            System.out.println("Invalid data was passed to the function matrix_from_array");
        }
        int nrow=arr.length/ncol;
        double[][] res = new double[nrow][];
        for(int i=0;i< nrow;i++) {
            double[] cur = new double[ncol];
            for(int j=0;j<ncol;j++) {
                cur[j]=arr[i+(j*nrow)];
            }
            res[i]=cur;
        }
        return res;
    }
}